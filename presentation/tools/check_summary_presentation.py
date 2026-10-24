#!/usr/bin/env python3
"""Check every structured-summary slide, diagnostic switch, and local artifact."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import hashlib
import json
from pathlib import Path
import threading
from urllib.parse import urljoin
from urllib.request import urlopen
from playwright.sync_api import sync_playwright


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    repo=Path(__file__).resolve().parents[2]
    class Handler(SimpleHTTPRequestHandler):
        def log_message(self,*args):pass
        def translate_path(self,path):
            if path.startswith('/offshell_reweighting/'):path=path[len('/offshell_reweighting'):]
            return super().translate_path(path)
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Handler,directory=str(repo)))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    origin=f'http://127.0.0.1:{server.server_port}';base=origin+'/offshell_reweighting/presentation/'
    report={'javascript_errors':[],'desktop_bounds':[],'local_links_checked':0,'diagnostic_views_checked':0}
    try:
        with sync_playwright() as pw:
            browser=pw.chromium.launch(headless=True,args=['--no-sandbox'])
            page=browser.new_page(viewport={'width':1440,'height':900},device_scale_factor=1)
            page.on('pageerror',lambda error:report['javascript_errors'].append(str(error)))
            page.goto(origin+'/offshell_reweighting/');page.wait_for_url('**/presentation/**')
            page.wait_for_function('window.Reveal && Reveal.isReady()')
            assert page.evaluate('window.SUMMARY_RESULTS.complete && window.FAMILY_RESULTS.complete')
            assert page.evaluate('FAMILY_RESULTS.models.every(m=>m.points===61) && FAMILY_RESULTS.grid[20]===1')
            ids=page.locator('.slides > section').evaluate_all('els=>els.map(e=>e.id)')
            assert len(ids)==33 and len(set(ids))==33
            assert ids[1:5]==['paper-result','same-sample','exact-closure','shape-systematics']
            assert ids[5:9]==['native-inputs','native-sr','native-channels','native-package']
            assert ids[9:12]==['normalizer-graph','normalizer-load','normalizer-create']
            assert ids[13:19]==['family-true-4096','family-data-like-4096','family-pruned-4096','family-true-8192','family-data-like-8192','family-pruned-8192']
            assert ids[20:24]==['pulls-4096','uncertainties-4096','pulls-8192','uncertainties-8192']
            assert page.locator('#resource-summary-table tbody tr').count()==6
            assert page.locator('[data-paper-slot] .paper-slot').count()==4
            assert page.locator('[data-paper-slot] img').count()==0
            for i,section_id in enumerate(ids):
                page.evaluate('i=>Reveal.slide(i)',i);page.wait_for_timeout(220)
                bounds=page.locator('#'+section_id).evaluate('''el=>{const rs=[...el.children].filter(e=>!e.hidden).map(e=>e.getBoundingClientRect());return {id:el.id,bottom:Math.max(...rs.map(r=>r.bottom)),right:Math.max(...rs.map(r=>r.right))};}''')
                report['desktop_bounds'].append(bounds)
                page.screenshot(path=str(args.output/f'slide_{i+1:02d}.png'))
                assert bounds['bottom']<870,bounds
                assert bounds['right']<1440,bounds
                assert 'NaN' not in page.locator('#'+section_id).inner_text()
            for size in (4096,8192):
                for kind in ('pulls','uncertainties'):
                    section_id=f'{kind}-{size}';page.evaluate('i=>Reveal.slide(i)',ids.index(section_id));page.wait_for_timeout(200)
                    buttons=page.locator(f'#{section_id} [data-diagnostic-view]')
                    for i in range(2):
                        button=buttons.nth(i);view=button.get_attribute('data-diagnostic-view');button.click()
                        suffix='pulls_'+view if kind=='pulls' else view
                        assert page.locator(f'#{section_id} [data-diagnostic-figure]').get_attribute('src').endswith(suffix+'.svg')
                        assert page.locator(f'#{section_id} [data-diagnostic-pdf]').get_attribute('href').endswith(suffix+'.pdf')
                        assert button.get_attribute('aria-pressed')=='true'
                        page.wait_for_function('sel=>{const e=document.querySelector(sel);return e.complete&&e.naturalWidth>0}',arg=f'#{section_id} [data-diagnostic-figure]')
                        page.screenshot(path=str(args.output/f'{section_id}_{view}.png'))
                        report['diagnostic_views_checked']+=1
                    buttons.nth(0).click()
            assert page.locator('img').evaluate_all('els=>els.every(e=>e.complete&&e.naturalWidth>0)')
            for path,digest in page.evaluate('Object.entries(SUMMARY_RESULTS.artifact_sha256)'):
                with urlopen(urljoin(base,path)) as response:assert hashlib.sha256(response.read()).hexdigest()==digest,path
            links=page.locator('[href],[src]').evaluate_all("els=>[...new Set(els.flatMap(e=>[e.getAttribute('href'),e.getAttribute('src')]).filter(Boolean))]")
            for link in links:
                if link.startswith(('#','https:','http:','data:','mailto:')):continue
                with urlopen(urljoin(base,link)) as response:assert response.status==200 and response.read(),link
                report['local_links_checked']+=1
            page.evaluate('Reveal.slide(0)');page.keyboard.press('ArrowRight');page.wait_for_timeout(250)
            assert page.evaluate('Reveal.getIndices().h')==1
            page.locator('#contents-button').click();page.locator('#contents-list a[href="#/cost-pruned-8192"]').click()
            assert page.evaluate('Reveal.getCurrentSlide().id')=='cost-pruned-8192'
            page.keyboard.press('Escape');assert page.evaluate('Reveal.isOverview()');page.keyboard.press('Escape')
            mobile=browser.new_page(viewport={'width':390,'height':844})
            mobile.on('pageerror',lambda error:report['javascript_errors'].append(str(error)))
            for section_id in ('paper-result','native-inputs','native-package','normalizer-create','family-true-4096','pulls-8192','cost-pruned-8192'):
                mobile.goto(base+'#/'+section_id);mobile.wait_for_url(base+'#'+section_id)
                assert mobile.locator('body').evaluate("el=>el.classList.contains('reading')")
                assert mobile.evaluate('document.documentElement.scrollWidth<=innerWidth'),'Mobile overflow'
                assert abs(mobile.locator('#'+section_id).bounding_box()['y']-48)<3
                mobile.screenshot(path=str(args.output/f'mobile_{section_id}.png'))
            mobile.locator('#contents-button').click();mobile.locator('#contents-list a[href="#uncertainties-4096"]').click()
            assert not mobile.locator('#contents').evaluate('el=>el.open')
            assert abs(mobile.locator('#uncertainties-4096').bounding_box()['y']-48)<3
            mobile.locator('#uncertainties-4096 [data-diagnostic-view="source_groups"]').click()
            assert mobile.evaluate('document.documentElement.scrollWidth<=innerWidth')
            # Verify configured paper slots without publishing a substitute paper figure.
            fixture=browser.new_page(viewport={'width':1440,'height':900})
            fixture.route('**/data/paper_figures.js',lambda route:route.fulfill(content_type='text/javascript',body='window.PAPER_FIGURES={pulls_4096:{src:"assets/summary_data_like_4096_pulls_nominal.svg",alt:"Layout test fixture",caption:"Test fixture"}};'))
            fixture.goto(base+'#/pulls-4096');fixture.wait_for_function('Reveal.isReady()')
            assert fixture.locator('[data-paper-slot="pulls_4096"] img').get_attribute('alt')=='Layout test fixture'
            assert fixture.locator('[data-paper-slot="pulls_4096"] figcaption').inner_text()=='Test fixture'
            page.goto(base+'history-20261023.html');page.wait_for_function('Reveal.isReady()')
            assert page.locator('.slides>section').count()==35
            assert page.evaluate('window.RESULTS.status === "complete" && window.FAMILY_RESULTS.complete')
            browser.close()
        assert not report['javascript_errors'],report['javascript_errors']
        report.update(status='passed',slides=33,navigation='passed',mobile='passed',paper_slots='passed',archive='passed')
    finally:
        server.shutdown();(args.output/'browser_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
