# auto_crispor.py
from pathlib import Path
from playwright.sync_api import sync_playwright, TimeoutError as PWTimeout

def run_auto_crispor(fasta_path: str, out_tsv: str, headless: bool = True):
    FASTA_PATH = Path(fasta_path)
    OUT_TSV = Path(out_tsv)

    if not FASTA_PATH.exists():
        raise FileNotFoundError(f"FASTA not found: {FASTA_PATH}")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        ctx = browser.new_context(accept_downloads=True)
        page = ctx.new_page()

        page.goto("https://crispor.gi.ucsc.edu/", wait_until="networkidle")

        # Find form frame (or main)
        frame = next((f for f in page.frames if f != page.main_frame and f.locator("textarea").first.count() > 0), page.main_frame)

        # Fill FASTA
        fasta_text = FASTA_PATH.read_text()
        try:
            ta = frame.locator("textarea[tabindex='1']").first
            if ta.count() == 0:
                ta = frame.locator("textarea").first
            ta.wait_for(timeout=10000)
            ta.fill(fasta_text)
        except PWTimeout:
            file_input = frame.locator("input[type='file']").first
            file_input.set_input_files(FASTA_PATH.as_posix())

        # Set genome (hidden select → JS)
        frame.locator("#genomeDropDown, select[name='org'], select[name='genome']").first.wait_for(state="attached", timeout=30000)
        frame.evaluate("""
          const sel = document.querySelector('#genomeDropDown, select[name="org"], select[name="genome"]');
          if (!sel) throw new Error('Genome <select> not found');
          let set = false;
          for (const opt of sel.options) {
            const t = (opt.text || '').toLowerCase();
            if (t.includes('homo sapiens') && (t.includes('grch38') || t.includes('hg38'))) {
              sel.value = opt.value;
              sel.dispatchEvent(new Event('change', { bubbles: true }));
              set = true;
              break;
            }
          }
          if (!set) throw new Error('Could not find Homo sapiens GRCh38/hg38 option');
        """)

        # Set PAM = NGG
        pam = None
        for sel in ["select[name='pam'][tabindex='3']", "select[name='pam']", "#pam"]:
            loc = frame.locator(sel)
            if loc.count() > 0:
                pam = loc
                break
        if pam:
            try:
                pam.select_option(value="NGG")
            except Exception:
                frame.evaluate("""
                  const s = document.querySelector('select[name="pam"], #pam');
                  let ok = false;
                  for (const o of s.options) {
                    const v = (o.value || '').toUpperCase();
                    const t = (o.text || '');
                    if (v === 'NGG' || t.includes('NGG')) {
                      s.value = o.value;
                      s.dispatchEvent(new Event('change', { bubbles: true }));
                      ok = true;
                      break;
                    }
                  }
                  if (!ok) throw new Error('NGG PAM not found');
                """)
        else:
            frame.evaluate("""
              const s = document.querySelector('select[name="pam"], #pam');
              if (!s) throw new Error('PAM <select> not found');
              let ok = false;
              for (const o of s.options) {
                const v = (o.value || '').toUpperCase();
                const t = (o.text || '');
                if (v === 'NGG' || t.includes('NGG')) {
                  s.value = o.value;
                  s.dispatchEvent(new Event('change', { bubbles: true }));
                  ok = true;
                  break;
                }
              }
              if (!ok) throw new Error('NGG PAM not found');
            """)

        # Submit
        submit_btn = frame.locator('input[type="submit"][name="submit"][tabindex="4"]')
        submit_btn.wait_for(timeout=10000)
        submit_btn.click()

        # Wait for and download TSV
        frame.wait_for_selector("a[href*='download=guides'][href*='format=tsv']", timeout=180000)
        with page.expect_download() as dl_info:
            frame.locator("a[href*='download=guides'][href*='format=tsv']").click()
        dl = dl_info.value
        dl.save_as(OUT_TSV.as_posix())

        ctx.close()
        browser.close()

# Allow CLI usage too
if __name__ == "__main__":
    # simple CLI: python auto_crispor.py <fasta> <out.tsv> [--show]
    import sys
    headless = True
    if len(sys.argv) >= 3:
        fasta = sys.argv[1]
        tsv = sys.argv[2]
        if len(sys.argv) >= 4 and sys.argv[3] == "--show":
            headless = False
        run_auto_crispor(fasta, tsv, headless=headless)
    else:
        # default for convenience
        run_auto_crispor("output_ENSG00000008086_3prime_sgRNAs_for_crispor.fasta",
                         "guides_hg38-unknownLoc.tsv",
                         headless=True)