#!/usr/bin/env python3
"""Gera a galeria HTML das figuras do PanNosoVax (imagens embutidas em base64).

Saída: results/report/gallery.html  (autocontida — publicável como Artifact)
"""
from __future__ import annotations
import base64, io, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from PIL import Image
from common import ROOT

FIGDIR = ROOT / "results/report/figures"
MAXW = 1500


def data_uri(name: str) -> str:
    im = Image.open(FIGDIR / name).convert("RGB")
    if im.width > MAXW:
        im = im.resize((MAXW, round(im.height * MAXW / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, "JPEG", quality=88)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def gif_uri(name: str) -> str:
    raw = (FIGDIR / name).read_bytes()
    return "data:image/gif;base64," + base64.b64encode(raw).decode()


FIGS = [
    ("F0_graphical_abstract.png", "Resumo", "Resumo gráfico completo", True,
     "Do genoma ao imunógeno numa figura: molécula direcionada aos três patógenos, fluxo do "
     "pipeline com as duas novidades, cobertura HLA e números-chave."),
    ("F1_funil_epitopos.png", "Figura 1", "Atrição de epitopos", False,
     "De dezenas de milhares de epitopos preditos a algumas dezenas selecionados — "
     "predição → conservação (≥95%) → segurança → cobertura."),
    ("F2_triagem_seguranca.png", "Figura 2", "Triagem negativa de segurança", False,
     "A camada de novidade: além de self humano e mimetismo, triagem contra o microbioma "
     "comensal. Em S. pneumoniae, 77% dos epitopos são barrados por homologia comensal."),
    ("F3_cobertura_populacional.png", "Figura 3", "Cobertura HLA — mundo vs Brasil", False,
     "O conjunto selecionado cobre ~96–98% da população, ponderado para frequências HLA "
     "brasileiras, com 4–6 epitopos por organismo."),
    ("F4_construto.png", "Figura 4", "Arquitetura do construto", False,
     "PanNosoVax_v1: 788 aa, adjuvante RS09, PADRE, 12 epitopos B + 15 MHC-II + 17 MHC-I, "
     "com linkers específicos por classe e sem cisteína livre."),
    ("F5_crossmatch_estrutural.png", "Figura 5", "Convergência estrutural", False,
     "Antígenos de superfície com a mesma dobra entre os três patógenos (TM-score até 0,89) — "
     "a dobra ABC-substrate-binding como alvo estrutural pan-patógeno."),
    ("F6_molecula.png", "Figura 6", "Modelo 3D do imunógeno", True,
     "Backbone do construto (modelo ESMFold) como tubo sombreado, colorido por bloco funcional, "
     "com o adjuvante RS09 em destaque. Versão de fundo claro disponível para artigo."),
]

PROV = ("F7_plddt.png", "Suplementar", "Confiança estrutural (pLDDT) — provisório", False,
        "pLDDT do modelo ESMFold single-sequence. Sai uniformemente baixo porque o ESMFold não "
        "dobra bem uma quimera de epitopos costurados — é limitação do método, não do construto. "
        "Avaliação definitiva aguarda o modelo ColabFold/AlphaFold (estágio 10).")

STATS = [("3", "patógenos"), ("302", "antígenos de superfície"), ("44", "epitopos no construto"),
         ("~98%", "cobertura HLA"), ("788 aa", "imunógeno · CAI 0,84")]


def card(fig, provisional=False):
    name, eyebrow, title, dark, cap = fig
    cls = "card" + (" card--dark" if dark else "") + (" card--prov" if provisional else "")
    banner = ('<p class="prov-tag">provisório · aguarda ColabFold</p>' if provisional else "")
    return f"""    <figure class="{cls}">
      <div class="card__head">
        <span class="eyebrow">{eyebrow}</span>
        <h3>{title}</h3>
      </div>
      {banner}
      <div class="card__img"><img loading="lazy" src="{data_uri(name)}" alt="{title}"></div>
      <figcaption>{cap}</figcaption>
    </figure>"""


def main():
    hero = data_uri("F8_pan_nosocomial.png")
    chips = "".join(f'<div class="chip"><b>{v}</b><span>{l}</span></div>' for v, l in STATS)
    cards = "\n".join(card(f) for f in FIGS)
    prov = card(PROV, provisional=True)
    gif_card = f"""    <figure class="card card--dark">
      <div class="card__head"><span class="eyebrow">Animação</span>
        <h3>Rotação 360° do imunógeno</h3></div>
      <div class="card__img"><img loading="lazy" src="{gif_uri('F6_molecula_360.gif')}" alt="molécula girando"></div>
      <figcaption>Rotação completa do modelo 3D — para apresentações e web.</figcaption>
    </figure>"""

    html = f"""<title>PanNosoVax — galeria de figuras</title>
<style>
:root{{
  --bg:#F4F6F9; --panel:#FFFFFF; --ink:#141A24; --sub:#5B6675; --line:#E2E7EE;
  --blue:#3B6EA5; --red:#D1495B; --teal:#2A9D8F; --amber:#D98A2B;
  --dark-bg:#0E1420; --shadow:0 1px 3px rgba(20,26,36,.08),0 8px 24px rgba(20,26,36,.06);
}}
@media (prefers-color-scheme:dark){{:root{{
  --bg:#0B1017; --panel:#141C28; --ink:#EAF0F7; --sub:#98A3B3; --line:#22303F;
  --shadow:0 1px 3px rgba(0,0,0,.4),0 10px 30px rgba(0,0,0,.35);
}}}}
:root[data-theme="dark"]{{
  --bg:#0B1017; --panel:#141C28; --ink:#EAF0F7; --sub:#98A3B3; --line:#22303F;
  --shadow:0 1px 3px rgba(0,0,0,.4),0 10px 30px rgba(0,0,0,.35);
}}
:root[data-theme="light"]{{
  --bg:#F4F6F9; --panel:#FFFFFF; --ink:#141A24; --sub:#5B6675; --line:#E2E7EE;
  --shadow:0 1px 3px rgba(20,26,36,.08),0 8px 24px rgba(20,26,36,.06);
}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);
  font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  line-height:1.55;-webkit-font-smoothing:antialiased}}
.wrap{{max-width:1040px;margin:0 auto;padding:clamp(20px,4vw,56px)}}
.mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}}
header{{margin-bottom:40px}}
.kicker{{font-family:ui-monospace,monospace;font-size:.72rem;letter-spacing:.22em;
  text-transform:uppercase;color:var(--teal);margin:0 0 14px}}
h1{{font-size:clamp(2.3rem,6vw,3.6rem);line-height:1.02;margin:0;font-weight:800;
  letter-spacing:-.02em;text-wrap:balance}}
h1 .em{{background:linear-gradient(100deg,var(--blue),var(--teal));-webkit-background-clip:text;
  background-clip:text;color:transparent}}
.lede{{color:var(--sub);font-size:1.12rem;max-width:60ch;margin:16px 0 0}}
.chips{{display:flex;flex-wrap:wrap;gap:10px;margin-top:26px}}
.chip{{display:flex;flex-direction:column;padding:10px 16px;border:1px solid var(--line);
  border-radius:12px;background:var(--panel);min-width:96px}}
.chip b{{font-size:1.3rem;font-weight:800;font-variant-numeric:tabular-nums}}
.chip span{{font-size:.74rem;color:var(--sub)}}
.hero{{margin:8px 0 48px;border-radius:18px;overflow:hidden;box-shadow:var(--shadow);
  background:var(--dark-bg)}}
.hero img{{width:100%;display:block}}
.grid{{display:flex;flex-direction:column;gap:34px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:16px;
  overflow:hidden;box-shadow:var(--shadow)}}
.card__head{{padding:20px 24px 4px}}
.eyebrow{{font-family:ui-monospace,monospace;font-size:.7rem;letter-spacing:.18em;
  text-transform:uppercase;color:var(--sub)}}
.card h3{{margin:6px 0 0;font-size:1.32rem;font-weight:750;letter-spacing:-.01em}}
.card__img{{padding:16px 24px 0}}
.card__img img{{width:100%;display:block;border-radius:8px}}
.card--dark .card__img{{background:var(--dark-bg);margin:16px 24px 0;padding:0;border-radius:10px}}
figcaption{{padding:14px 24px 22px;color:var(--sub);font-size:.95rem;max-width:66ch}}
.rule{{display:flex;align-items:center;gap:14px;margin:52px 0 26px;color:var(--sub);
  font-family:ui-monospace,monospace;font-size:.72rem;letter-spacing:.16em;text-transform:uppercase}}
.rule::before,.rule::after{{content:"";height:1px;background:var(--line);flex:1}}
.card--prov{{border-style:dashed;opacity:.94}}
.prov-tag{{margin:0;padding:8px 24px;color:var(--amber);font-size:.78rem;
  font-family:ui-monospace,monospace;letter-spacing:.05em}}
footer{{margin-top:56px;padding-top:24px;border-top:1px solid var(--line);color:var(--sub);
  font-size:.86rem}}
footer .mono{{color:var(--ink)}}
a{{color:var(--blue)}}
</style>

<div class="wrap">
  <header>
    <p class="kicker">Vacinologia reversa · desenho in silico</p>
    <h1><span class="em">PanNosoVax</span></h1>
    <p class="lede">Vacina multi-epitopo contra pneumonia nosocomial, cobrindo três patógenos
      prioritários da OMS — <i>Klebsiella pneumoniae</i>, <i>Acinetobacter baumannii</i> e
      <i>Streptococcus pneumoniae</i> — com um único imunógeno quimérico.</p>
    <div class="chips">{chips}</div>
  </header>

  <div class="hero"><img src="{hero}" alt="Resumo gráfico do PanNosoVax"></div>

  <div class="grid">
{cards}
{gif_card}
  </div>

  <div class="rule">material suplementar</div>
  <div class="grid">
{prov}
  </div>

  <footer>
    Figuras geradas a partir dos dados reais do pipeline · paleta validada para daltonismo ·
    PNG 300 dpi + PDF + SVG em <span class="mono">results/report/figures/</span>.
    Cada figura segue a regra do projeto: <b>figura sem dado não é desenhada</b>.
  </footer>
</div>
"""
    out = ROOT / "results/report/gallery.html"
    out.write_text(html)
    print(f"galeria escrita: {out} ({len(html)//1024} KB)")


if __name__ == "__main__":
    main()
