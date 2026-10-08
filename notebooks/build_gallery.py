"""Build a local, self-contained navigation page from actually exported figures."""
from pathlib import Path
from html import escape

directory = Path(__file__).resolve().parent
topics = [
    ("01_tracking_concepts", "3-D tracking and boundaries", "MATLAB-verified periodic boundaries, splits, and merges; two different 3-D event definitions."),
    ("02_real_3d", "Real three-dimensional data", "400×160×731: occurrence frequency, event statistics, and real event evolution."),
    ("03_regular_4d", "Regular four-dimensional connectivity", "Periodic longitude and depth/time diagonals; exact event membership comparisons."),
    ("04_llc_seams", "Actual LLC90 topology", "Locations, orientation, and exhaustive connectivity checks for 2,312 seams."),
    ("05_auxiliary_spn", "Smoothing and normalization", "Spatial/temporal sampling of real events and numeric residuals against MATLAB."),
    ("06_validation_performance", "Validation and performance", "Recomputed MATLAB fixtures, recorded 3-D timings, and current 4-D measurements."),
    ("07_tasman_sea_spn", "The full Tasman Sea normalization example", "230 days of SST anomaly, both norm_flag modes, five normalized times, and three complete GIF animations."),
]
cards = []
for stem, title, description in topics:
    images = sorted((directory / "figures").glob(stem[:2] + "_*.png"))
    figures = "".join(
        f'<figure><a href="figures/{escape(p.name)}"><img loading="lazy" src="figures/{escape(p.name)}" alt="{escape(p.stem)}"></a>'
        f'<figcaption>{escape(p.stem)} · <a href="figures/{escape(p.with_suffix(".pdf").name)}">PDF</a></figcaption></figure>'
        for p in images)
    animations = sorted((directory / "animations").glob(stem[:2] + "_*.gif"))
    animation_links = "".join(f'<a href="animations/{escape(p.name)}">{escape(p.stem)} · GIF</a>' for p in animations)
    cards.append(f'<section><div class="topic"><span>{stem[:2]}</span><h2>{title}</h2></div>'
                 f'<p>{description}</p><nav><a href="html/{stem}.html">View executed notebook</a>'
                 f'<a href="{stem}.ipynb">Jupyter file</a>{animation_links}</nav><div class="figures">{figures}</div></section>')
page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>MHWTracking · Notebooks & Figures</title><style>
*{box-sizing:border-box}body{margin:0;background:#f3f6f9;color:#162635;font-family:system-ui,-apple-system,sans-serif;line-height:1.65}
header{background:#102e45;color:white;padding:48px max(24px,calc((100vw - 1250px)/2)) 36px}header p{max-width:850px;color:#cddfe9}
h1{font-size:34px;letter-spacing:-.5px;margin:4px 0}main{max-width:1300px;margin:auto;padding:28px 24px}section{background:white;margin-bottom:28px;padding:24px;border:1px solid #dde5eb;border-radius:14px}
.topic{display:flex;align-items:center;gap:14px}.topic span{background:#e4f3f2;color:#066f71;font-size:22px;font-weight:700;padding:8px 14px;border-radius:8px}h2{margin:0;font-size:24px}nav{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}a{color:#056977}nav a{background:#edf6f6;padding:8px 16px;text-decoration:none;border-radius:7px}.figures{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(440px,100%),1fr));gap:18px}figure{margin:0;border:1px solid #e4e9ed;border-radius:8px;overflow:hidden}img{width:100%;display:block;background:white}figcaption{font-size:12px;padding:9px 12px;background:#f8fafb;color:#52616d}footer{padding:10px 0 30px;color:#52616d;font-size:14px}.tag{font-size:13px;letter-spacing:2px;text-transform:uppercase;color:#8bd3ce}
</style><header><div class="tag">Julia · MATLAB validation · ocean events</div><h1>MHWTracking · Notebooks & Figures</h1><p>Six runnable Julia notebooks with figures from actual execution and MATLAB validation assertions. Each notebook identifies real 3-D data, synthetic 4-D cases, actual LLC topology, and recorded performance measurements.</p></header><main>'''
page = page.replace("Six runnable Julia notebooks", "Seven runnable Julia notebooks")
page += "".join(cards)
page += '<footer>See README.md for the notebook environment and reproduction commands. Each HTML notebook embeds code outputs and figures. Scientific definitions and validation scope are documented in the project docs.</footer></main></html>'
(directory / "gallery.html").write_text(page)
print(f"Gallery: {len(topics)} topics, {sum(len(list((directory / 'figures').glob(stem[:2] + '_*.png'))) for stem,_,_ in topics)} PNG figures")
