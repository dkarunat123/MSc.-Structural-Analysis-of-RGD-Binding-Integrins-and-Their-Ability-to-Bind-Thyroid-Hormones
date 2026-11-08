from pathlib import Path
import csv,json
base=Path(__file__).resolve().parents[1] / 'docs'
with (base/'data/receptor_affinities.csv').open(encoding='utf8',newline='') as f: rec=list(csv.DictReader(f))
for x in rec:
 x['neg_log10_kd']=float(x['neg_log10_kd']);x['kd_m']=float(x['kd_m'])
with (base/'data/species_medians.csv').open(encoding='utf8',newline='') as f: med=list(csv.DictReader(f))
for x in med:
 x['median_neg_log10_kd']=float(x['median_neg_log10_kd']);x['n_receptors']=int(x['n_receptors'])
with (base/'data/cluster_orders.json').open(encoding='utf8') as f: clusters=json.load(f)
with (base/'data/species_rank_correlations.csv').open(encoding='utf8',newline='') as f: spcorr=list(csv.DictReader(f))
for x in spcorr:x['spearman_rho']=float(x['spearman_rho'])
with (base/'assets/data.js').open('w',encoding='utf8') as f:
 f.write('/* Data are bundled for offline use. See docs/data/validation_report.json. */\n')
 f.write('window.THESIS_DATA = '+json.dumps({'receptors':rec,'medians':med,'clusters':clusters,'speciesCorrelations':spcorr},ensure_ascii=False,separators=(',',':'))+';\n')
