(() => {
  'use strict';
  // Version: sea-urchin-cell-numbers-v2 (2026-11-08). All 24 rows show their values.
  const D = window.THESIS_DATA;
  const SPECIES = ['H. sapiens', 'S. purpuratus', 'N. vectensis'];
  const LIGANDS = ['Cilengitide','T3','T4','TRIAC','TETRAC','sT3','sT4','rT3','T4AM'];
  const PALETTE_AFF = ['#313695','#74ADD1','#FFFFBF','#FDAE61','#A50026'];
  const PALETTE_RHO = ['#313695','#74ADD1','#FFFFBF','#F46D43','#A50026'];
  // Thesis-wide structure/model identity colours; heatmap palettes above encode numeric scores.
  const MODEL_COLORS = {'H. sapiens':'#396ea5','S. purpuratus':'#a641a6','N. vectensis':'#6f990c'};
  const MODEL_LIGHT = {'H. sapiens':'#8bc1f0','S. purpuratus':'#ec9de5','N. vectensis':'#c9df7d'};
  const $ = q => document.querySelector(q);
  const all = q => [...document.querySelectorAll(q)];
  const NS='http://www.w3.org/2000/svg';
  const recIndex=new Map(D.receptors.map(x=>[`${x.species}|${x.receptor}|${x.ligand}`,x]));
  const medianIndex=new Map(D.medians.map(x=>[`${x.species}|${x.ligand}`,x]));
  let figure='14', currentChart=null, selectedId=null;
  let rhoLookup={};
  for(const sp of SPECIES){
    const recIds=D.clusters[sp].figure14.row_ids;
    const ranks=LIGANDS.map(l=>recIds.map(r=>recIndex.get(`${sp}|${r}|${l}`).neg_log10_kd));
    const rankColumns=ranks.map(rankValues);
    rhoLookup[sp]=LIGANDS.map((_,i)=>LIGANDS.map((_,j)=>pearson(rankColumns[i],rankColumns[j])));
  }
  function rankValues(v){const sorted=[...v].map((n,i)=>({n,i})).sort((a,b)=>a.n-b.n);let output=Array(v.length);for(let i=0;i<sorted.length;){let j=i+1;while(j<sorted.length&&sorted[j].n===sorted[i].n)j++;const rank=(i+1+j)/2;for(let k=i;k<j;k++)output[sorted[k].i]=rank;i=j;}return output;}
  function pearson(a,b){const ma=mean(a),mb=mean(b);let n=0,da=0,db=0;for(let i=0;i<a.length;i++){const x=a[i]-ma,y=b[i]-mb;n+=x*y;da+=x*x;db+=y*y;}return n/Math.sqrt(da*db);}
  function mean(a){return a.reduce((s,n)=>s+n,0)/a.length;}
  function median(a){const x=[...a].sort((a,b)=>a-b),n=x.length;return n%2?x[(n-1)/2]:(x[n/2-1]+x[n/2])/2;}
  function elem(name,attrs={},text){const x=document.createElementNS(NS,name);for(const [key,val] of Object.entries(attrs))x.setAttribute(key,String(val));if(text!==undefined)x.textContent=String(text);return x;}
  function cssColor(hex){const s=hex.replace('#','');return [0,2,4].map(i=>parseInt(s.slice(i,i+2),16));}
  function interpolate(value,min,max,palette){const t=Math.max(0,Math.min(0.999999,(value-min)/((max-min)||1)))*(palette.length-1),i=Math.floor(t),f=t-i;const a=cssColor(palette[i]),b=cssColor(palette[i+1]);return '#'+a.map((v,k)=>Math.round(v*(1-f)+b[k]*f).toString(16).padStart(2,'0')).join('');}
  function inkOn(bg){const [r,g,b]=cssColor(bg);return (0.2126*r+0.7152*g+0.0722*b)<145?'#ffffff':'#132235';}
  function chartCellId(row,col){return `${figure}|${row}|${col}`;}
  function hint(e,markup){const tip=$('#tooltip');tip.innerHTML=markup;tip.style.display='block';const x=Math.min(innerWidth-345,e.clientX+17),y=Math.min(innerHeight-125,e.clientY+13);tip.style.left=Math.max(10,x)+'px';tip.style.top=Math.max(10,y)+'px';}
  function hideHint(){$('#tooltip').style.display='none';}
  function escape(s){return String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
  function saveFile(name,blob){const a=document.createElement('a');const url=URL.createObjectURL(blob);a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1200);}
  function clusterSegments(z,items,xPositions,yBottom,yTop){
    const n=items.length,heights=z.map(r=>r[2]),maxH=Math.max(...heights)||1;
    const nodes=Array(n+z.length);
    for(let i=0;i<n;i++)nodes[i]={x:xPositions.get(i),y:yBottom};
    const lines=[];
    for(let i=0;i<z.length;i++){
      const [a,b,h]=z[i],v1=nodes[a],v2=nodes[b];
      const y=yBottom-(h/maxH)*(yBottom-yTop);
      if(!v1||!v2||!Number.isFinite(v1.x)||!Number.isFinite(v2.x))continue;
      lines.push([v1.x,v1.y,v1.x,y],[v2.x,v2.y,v2.x,y],[v1.x,y,v2.x,y]);
      nodes[n+i]={x:(v1.x+v2.x)/2,y};
    }
    return lines;
  }
  function getSpec(){
    const sp=$('#species').value, ordering=$('#ordering').value,filter=$('#search').value.trim().toLowerCase();
    let rows,cols=[],values,min,max,tree,kind,metadata={};
    if(figure==='14'){
      const ids=D.clusters[sp].figure14.row_ids;
      const z=D.clusters[sp].figure14;
      let orderedRow=ids.slice(), orderedLig=LIGANDS.slice();
      if(ordering==='cluster'){
        orderedRow=z.row_order.map(i=>ids[i]);orderedLig=z.ligand_order.map(i=>LIGANDS[i]);
      }else if(ordering==='alpha'){
        orderedRow.sort((a,b)=>recIndex.get(`${sp}|${a}|T3`).receptor_label.localeCompare(recIndex.get(`${sp}|${b}|T3`).receptor_label));orderedLig.sort((a,b)=>a.localeCompare(b));
      }else{
        orderedRow.sort((a,b)=>mean(LIGANDS.map(l=>recIndex.get(`${sp}|${b}|${l}`).neg_log10_kd))-mean(LIGANDS.map(l=>recIndex.get(`${sp}|${a}|${l}`).neg_log10_kd)));
        orderedLig.sort((a,b)=>median(ids.map(r=>recIndex.get(`${sp}|${r}|${b}`).neg_log10_kd))-median(ids.map(r=>recIndex.get(`${sp}|${r}|${a}`).neg_log10_kd)));
      }
      orderedLig=orderedLig.filter(l=>l.toLowerCase().includes(filter));
      rows=orderedRow.map(r=>({key:r,label:recIndex.get(`${sp}|${r}|T3`).receptor_label,sub:r}));
      cols=orderedLig.map(l=>({key:l,label:l}));
      values=rows.map(r=>cols.map(c=>recIndex.get(`${sp}|${r.key}|${c.key}`).neg_log10_kd));
      let allValues=D.receptors.filter(x=>$('#sharedScale').checked||x.species===sp).map(x=>x.neg_log10_kd);
      min=Math.min(...allValues);max=Math.max(...allValues);
      tree=ordering==='cluster'&&!filter?z.ligand_linkage:null;
      kind='affinity';metadata={sp,originalRowCount:ids.length};
    }else if(figure==='15'){
      const z=D.clusters[sp].figure15;
      let items=LIGANDS.map((_,i)=>i);
      if(ordering==='cluster')items=z.ligand_order.slice();
      else if(ordering==='alpha')items.sort((a,b)=>LIGANDS[a].localeCompare(LIGANDS[b]));
      else items.sort((a,b)=>mean(rhoLookup[sp][b].filter((_,j)=>j!==b))-mean(rhoLookup[sp][a].filter((_,j)=>j!==a)));
      rows=items.map(i=>({key:LIGANDS[i],label:LIGANDS[i]}));cols=items.map(i=>({key:LIGANDS[i],label:LIGANDS[i]}));
      values=rows.map(r=>cols.map(c=>rhoLookup[sp][LIGANDS.indexOf(r.key)][LIGANDS.indexOf(c.key)]));
      min=-1;max=1;tree=ordering==='cluster'?z.ligand_linkage:null;kind='rho';metadata={sp};
    }else{
      const z=D.clusters.figure16;
      let spOrder=SPECIES.slice(),ligOrder=LIGANDS.slice();
      if(ordering==='cluster'){spOrder=z.species_order.map(i=>SPECIES[i]);ligOrder=z.ligand_order.map(i=>LIGANDS[i]);}
      else if(ordering==='alpha'){spOrder.sort();ligOrder.sort((a,b)=>a.localeCompare(b));}
      else ligOrder.sort((a,b)=>mean(SPECIES.map(sp=>medianIndex.get(`${sp}|${b}`).median_neg_log10_kd))-mean(SPECIES.map(sp=>medianIndex.get(`${sp}|${a}`).median_neg_log10_kd)));
      ligOrder=ligOrder.filter(l=>l.toLowerCase().includes(filter));
      rows=spOrder.map(sp=>({key:sp,label:sp}));cols=ligOrder.map(l=>({key:l,label:l}));
      values=rows.map(r=>cols.map(c=>medianIndex.get(`${r.key}|${c.key}`).median_neg_log10_kd));
      min=Math.min(...D.medians.map(x=>x.median_neg_log10_kd));max=Math.max(...D.medians.map(x=>x.median_neg_log10_kd));
      tree=ordering==='cluster'&&!filter?z.ligand_linkage:null;kind='affinity';metadata={};
    }
    return {rows,cols,values,min,max,tree,kind,metadata};
  }
  function updateTitle(spec){
    const sp=spec.metadata.sp||'';const nice=({'H. sapiens':'Human','S. purpuratus':'Sea urchin','N. vectensis':'Sea anemone'})[sp];
    $('#figureId').textContent=`Figure ${figure} · interactive exploration`;
    if(figure==='14'){
      $('#figureTitle').textContent=`${nice} integrin binding-affinity profiles`;
      $('#figureDescription').textContent='Predicted -log₁₀(KD) values for nine ligands across individual integrin models. Default clustering uses Euclidean distance and complete linkage; each species has an independent colour scale, as in the thesis.';
      $('#method').innerHTML=`<b>${spec.rows.length} modeled receptors</b> × ${spec.cols.length} ligands · Euclidean / complete linkage`;
    }else if(figure==='15'){
      $('#figureTitle').textContent=`${nice} ligand-ligand similarities`;
      $('#figureDescription').textContent='Pairwise Spearman ρ values compare how two ligands rank receptor models within the selected species. These are not correlations of measurements from independent animals.';
      $('#method').innerHTML=`<b>9 × 9 correlations</b> · cluster distance 1 - ρ · complete linkage`;
    }else{
      $('#figureTitle').textContent='Cross-species median ligand-affinity profiles';
      $('#figureDescription').textContent='Each cell is the median predicted -log₁₀(KD) across the complete receptor repertoire modeled for that species (human 8, sea urchin 24, sea anemone 16).';
      $('#method').innerHTML='<b>3 species × 9 ligands</b> · unstandardized medians · Euclidean / complete linkage';
    }
    // Model identity is separate from the continuous affinity/correlation colour scale.
    const selectedSpecies = sp || 'H. sapiens';
    document.documentElement.style.setProperty('--active-model', figure==='16'?'#53627a':MODEL_COLORS[selectedSpecies]);
    document.documentElement.style.setProperty('--active-model-deep', figure==='16'?'#53627a':MODEL_COLORS[selectedSpecies]);
    $('#modelContext').textContent = figure==='16' ? 'All three species · AlphaFold models' : (nice + ' · AlphaFold model');
    $('#legendMin').textContent=spec.min.toFixed(2);$('#legendMax').textContent=spec.max.toFixed(2);
    $('#gradient').classList.toggle('rho',spec.kind==='rho');
    $('#legendDescription').textContent=spec.kind==='rho'?'Spearman ρ · -1 opposite ranks; +1 matching ranks':'Predicted -log₁₀(KD) · higher = stronger predicted binding';
    $('#speciesControl').style.display=figure==='16'?'none':'';
    $('#scaleControl').style.display=figure==='14'?'':'none';
    $('#searchControl').style.display=figure==='15'?'none':'';
    $('#ordering').options[2].textContent=figure==='15'?'Mean correlation: highest first':'Affinity: strongest first';
  }
  function drawTree(svg,spec,left,cellW){
    if(!spec.tree||spec.cols.length!==9)return;
    const origIds=spec.cols.map(c=>LIGANDS.indexOf(c.key));
    const positions=new Map(origIds.map((id,i)=>[id,left+i*cellW+cellW/2]));
    const lines=clusterSegments(spec.tree,LIGANDS,positions,90,17);
    const g=elem('g',{'aria-hidden':true});
    for(const [x1,y1,x2,y2] of lines){g.appendChild(elem('line',{x1,y1,x2,y2,stroke:figure==='16'?'#a6b1c1':MODEL_COLORS[spec.metadata.sp],'stroke-width':1.5,opacity:.85}));}svg.appendChild(g);
  }
  function selectCell(spec,row,col,value,tile){
    selectedId=chartCellId(row.key,col.key);
    const title=$('#cellTitle'),desc=$('#cellDescription'),score=$('#cellScore'),unit=$('#cellUnit');
    if(figure==='14'){
      const ent=recIndex.get(`${spec.metadata.sp}|${row.key}|${col.key}`);
      title.textContent=`${ent.receptor_label} × ${col.key} · ${spec.metadata.sp}`;
      desc.textContent=`Computational model ${ent.receptor}; predicted KD = ${ent.kd_m.toExponential(3)} M. Subunits ${ent.alpha_subunit} + ${ent.beta_subunit}; docking prediction, not a measured equilibrium constant.`;
      score.textContent=value.toFixed(3);unit.textContent='predicted -log₁₀(KD)';
    }else if(figure==='15'){
      title.textContent=`${row.key} ↔ ${col.key} · ${spec.metadata.sp}`;
      desc.textContent=`Spearman rank correlation across ${D.clusters[spec.metadata.sp].figure14.row_ids.length} modelled receptors. A positive value indicates similar receptor rank ordering for these ligands; it does not measure agreement of absolute affinity magnitudes.`;
      score.textContent=value.toFixed(3);unit.textContent='Spearman ρ';
    }else{
      const ent=medianIndex.get(`${row.key}|${col.key}`);
      title.textContent=`${col.key} · ${row.key}`;
      desc.textContent=`Median predicted affinity across ${ent.n_receptors} different computational heterodimer models in this species. The median is taken directly from the -log₁₀(KD) scores, without re-transforming raw KD.`;
      score.textContent=value.toFixed(3);unit.textContent='median -log₁₀(KD)';
    }
    all('.tile').forEach(x=>x.setAttribute('stroke-width','2'));
    all('.tile').forEach(x=>x.classList.remove('selectedtile'));
    if(tile){tile.classList.add('selectedtile');tile.setAttribute('stroke-width','3.5');}
  }
  function render(){
    const spec=getSpec();currentChart=spec;updateTitle(spec);
    const target=$('#chart');target.replaceChildren();
    if(!spec.cols.length){target.textContent='No ligand matches your search. Clear the filter to display the heatmap.';return;}
    const isCorr=figure==='15';
    const cw=isCorr?72:(figure==='16'?84:72),ch=isCorr?47:(figure==='16'?71:(spec.rows.length>18?28:36));
    const left=figure==='14'?163:178,top=145,right=30,bottom=30;
    const W=left+cw*spec.cols.length+right,H=top+ch*spec.rows.length+bottom;
    const svg=elem('svg',{xmlns:NS,viewBox:`0 0 ${W} ${H}`,width:W,height:H,role:'img','aria-label':`Figure ${figure} interactive heatmap with ${spec.rows.length} rows and ${spec.cols.length} columns`});
    const background=elem('rect',{x:left,y:top,width:cw*spec.cols.length,height:ch*spec.rows.length,fill:'#f3f5f9',rx:5});svg.appendChild(background);
    if(spec.tree)drawTree(svg,spec,left,cw);
    for(let j=0;j<spec.cols.length;j++){
      const x=left+j*cw+cw/2;
      svg.appendChild(elem('text',{x,y:117,'text-anchor':'middle',fill:'#405069','font-size':12,'font-weight':750},spec.cols[j].label));
    }
    for(let i=0;i<spec.rows.length;i++){
      const y=top+i*ch;
      const r=spec.rows[i];
      const modelSpecies = figure==='16'?r.key:spec.metadata.sp;
      const modelTint = MODEL_COLORS[modelSpecies] || '#526178';
      // Species-coloured rail identifies the model; the heatmap cell fill remains quantitative.
      svg.appendChild(elem('rect',{x:left-7,y:y+2,width:4,height:ch-4,rx:2,fill:MODEL_COLORS[modelSpecies]||'#60728c'}));
      svg.appendChild(elem('text',{x:left-14,y:y+ch/2+4,'text-anchor':'end',fill:modelTint,'font-size':figure==='14'?11:12,'font-weight':650},r.label));
      for(let j=0;j<spec.cols.length;j++){
        const c=spec.cols[j],v=spec.values[i][j];
        const x=left+j*cw;const base=interpolate(v,spec.min,spec.max,isCorr?PALETTE_RHO:PALETTE_AFF);
        const rect=elem('rect',{x:x+1,y:y+1,width:cw-2,height:ch-2,rx:3,fill:base,stroke:'#f1f3f7','stroke-width':2,tabindex:0,role:'button','aria-label':`${r.label}, ${c.label}: ${v.toFixed(4)}`});
        rect.classList.add('tile');
        if(selectedId===chartCellId(r.key,c.key)){rect.classList.add('selectedtile');rect.setAttribute('stroke-width','3.5');}
        rect.addEventListener('mouseenter',e=>hint(e,`<strong>${escape(r.label)} × ${escape(c.label)}</strong><br>${v.toFixed(4)} ${isCorr?'Spearman ρ':'predicted -log₁₀(KD)'}<br><span style="color:#62728a">Click for interpretation</span>`));
        rect.addEventListener('mousemove',e=>hint(e,`<strong>${escape(r.label)} × ${escape(c.label)}</strong><br>${v.toFixed(4)} ${isCorr?'Spearman ρ':'predicted -log₁₀(KD)'}<br><span style="color:#62728a">Click for interpretation</span>`));
        rect.addEventListener('mouseleave',hideHint);
        rect.addEventListener('click',()=>selectCell(spec,r,c,v,rect));
        rect.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectCell(spec,r,c,v,rect);}});
        svg.appendChild(rect);
        // Display numerical labels for every species, including the 24-row sea urchin heatmap.
        // Smaller cells use a slightly smaller font rather than hiding their values.
        svg.appendChild(elem('text',{x:x+cw/2,y:y+ch/2+4,'text-anchor':'middle',fill:inkOn(base),'font-size':ch<=30?9:ch<=38?10:12,'font-weight':780,'pointer-events':'none'},v.toFixed(2)));
      }
    }
    target.appendChild(svg);
    // Keep the original selected cell when switching ordering, if it is still visible.
    if(!selectedId || !all('.selectedtile').length){
      const r=spec.rows[0],c=spec.cols[0],v=spec.values[0][0];
      selectCell(spec,r,c,v,svg.querySelector('.tile'));
    }
  }
  function tableCSV(spec){const label=figure==='14'?'receptor':figure==='15'?'ligand':'species';const output=[[label,...spec.cols.map(x=>x.key)]];for(let i=0;i<spec.rows.length;i++)output.push([spec.rows[i].key,...spec.values[i].map(n=>n.toPrecision(15))]);return output.map(row=>row.map(x=>'"'+String(x).replace(/"/g,'""')+'"').join(',')).join('\r\n')+'\r\n';}
  function switchFigure(num){figure=num;selectedId=null;all('.figBtn').forEach(b=>{const yes=b.dataset.figure===num;b.classList.toggle('selected',yes);b.setAttribute('aria-pressed',String(yes));});$('#search').value='';render();}
  all('.figBtn').forEach(b=>b.addEventListener('click',()=>switchFigure(b.dataset.figure)));
  for(const id of ['species','ordering','sharedScale'])$('#'+id).addEventListener('change',()=>{selectedId=null;render();});
  $('#search').addEventListener('input',()=>{selectedId=null;render();});
  $('#exportCsv').addEventListener('click',()=>{if(currentChart)saveFile(`figure${figure}_${figure==='16'?'cross_species':$('#species').value.replace(/[^a-z]/gi,'_')}.csv`,new Blob([tableCSV(currentChart)],{type:'text/csv;charset=utf-8'}));});
  $('#exportSvg').addEventListener('click',()=>{const svg=$('#chart svg');if(!svg)return;const out='<?xml version="1.0" encoding="UTF-8"?>\n'+new XMLSerializer().serializeToString(svg);saveFile(`figure${figure}_interactive_export.svg`,new Blob([out],{type:'image/svg+xml;charset=utf-8'}));});
  render();
})();
