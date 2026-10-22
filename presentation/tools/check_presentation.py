#!/usr/bin/env python3
"""Check the rendered slideshow, mobile reading mode, and local download links."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import threading
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[2]

    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def translate_path(self, path):
            if path.startswith("/offshell_reweighting/"):
                path = path[len("/offshell_reweighting"):]
            return super().translate_path(path)

    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(Handler, directory=str(repo)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    origin = f"http://127.0.0.1:{server.server_port}"
    page_url = origin+"/offshell_reweighting/presentation/"
    report = {"javascript_errors": [], "desktop_bounds": [], "local_links_checked": 0}
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(headless=True, args=["--no-sandbox"])
            page = browser.new_page(viewport={"width": 1440, "height": 900}, device_scale_factor=1)
            page.on("pageerror", lambda error: report["javascript_errors"].append(str(error)))
            page.goto(origin+"/offshell_reweighting/")
            page.wait_for_url("**/presentation/**")
            page.wait_for_function("window.Reveal && Reveal.isReady()")
            assert page.evaluate("window.RESULTS.status") == "complete"
            ids = page.locator(".slides > section").evaluate_all("els => els.map(e => e.id)")
            for i, id_ in enumerate(ids):
                page.evaluate("i => Reveal.slide(i)", i)
                page.wait_for_timeout(300)
                bounds = page.locator(f"#{id_}").evaluate("""el => {
                  const rects = [...el.children].filter(e => !e.hidden).map(e => e.getBoundingClientRect());
                  return {id: el.id, bottom: Math.max(...rects.map(r=>r.bottom)), right: Math.max(...rects.map(r=>r.right))};
                }""")
                report["desktop_bounds"].append(bounds)
                page.screenshot(path=str(args.output/f"slide_{i+1:02d}.png"))
                assert bounds["bottom"] < 870, bounds
                assert bounds["right"] < 1440, bounds
            assert page.locator("img").evaluate_all("els => els.every(e => e.complete && e.naturalWidth > 0)")
            links = page.locator("[href],[src]").evaluate_all("els => [...new Set(els.flatMap(e=>[e.getAttribute('href'),e.getAttribute('src')]).filter(Boolean))]")
            for link in links:
                if link.startswith(("#", "http:", "https:", "data:", "mailto:")):
                    continue
                with urlopen(urljoin(page_url, link)) as response:
                    assert response.status == 200, link
                    assert response.read(), link
                report["local_links_checked"] += 1
            page.evaluate("Reveal.slide(0)")
            page.keyboard.press("ArrowRight")
            page.wait_for_timeout(250)
            assert page.evaluate("Reveal.getIndices().h") == 1
            page.locator("#contents-button").click()
            page.locator('#contents-list a[href="#/anchor-correction"]').click()
            assert page.evaluate("Reveal.getCurrentSlide().id") == "anchor-correction"
            page.keyboard.press("Escape")
            assert page.evaluate("Reveal.isOverview()")
            page.keyboard.press("Escape")
            mobile = browser.new_page(viewport={"width": 390, "height": 844}, device_scale_factor=1)
            mobile.on("pageerror", lambda error: report["javascript_errors"].append(str(error)))
            mobile.goto(page_url)
            assert mobile.locator("body").evaluate("el => el.classList.contains('reading')")
            assert mobile.evaluate("document.documentElement.scrollWidth <= innerWidth"), "Horizontal mobile overflow"
            mobile.locator("#contents-button").click()
            mobile.locator('#contents-list a[href="#scan-8192"]').click()
            mobile.wait_for_timeout(100)
            mobile.screenshot(path=str(args.output/"mobile_scan.png"))
            assert not mobile.locator("#contents").evaluate("el => el.open")
            assert abs(mobile.locator("#scan-8192").bounding_box()["y"]-48) < 3
            browser.close()
            assert not report["javascript_errors"], report["javascript_errors"]
            report.update(status="passed", slides=len(ids), navigation="passed", mobile_navigation="passed", project_subpath_redirect="passed")
    finally:
        server.shutdown()
        (args.output/"browser_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
