/* Additional offline report views; all inference results come from saved tables. */
(() => {
"use strict";
const C=window.proteomicsPreview,D=C.data,E=D.extra;
const {el,esc,finite,fmt,contrast,stats,statFor,proteinMap,sampleMap,groupColor,directionColor,config,axis,baseLayout,smallStat}=C;
const S={ora:null,gsea:null,gsva:null,ppi:null};
const diverging=[[0,"#2166ac"],[.5,"#f7f7f7"],[1,"#b35806"]];
const sequential=[[0,"#dceaf4"],[1,"#006ba4"]];
function fitLabels(layout,id,names){
  const narrow=el(id).clientWidth>0&&el(id).clientWidth<600;
  layout.margin.l=narrow?140:230;
  const limit=narrow?21:36;
  layout.yaxis.tickmode="array";layout.yaxis.tickvals=names;
  layout.yaxis.ticktext=names.map(n=>n.length>limit?n.slice(0,limit-1)+"…":n);
  layout.yaxis.tickfont={size:11};layout.yaxis.automargin=false;
}
const toGenes=s=>s?String(s).split(";").filter(Boolean):[];
const coName=()=>C.getState().contrast;
function variance(v){const a=v.filter(finite);if(a.length<2)return -1;const m=a.reduce((s,x)=>s+x,0)/a.length;return a.reduce((s,x)=>s+(x-m)**2,0)/(a.length-1);}
function geneProteins(genes){const set=new Set(genes);return D.proteins.filter(p=>set.has(p.gene));}
function proteinTable(genes){
  const pp=geneProteins(genes);
  return "<table><thead><tr><th>基因</th><th>蛋白 ID</th><th>log2FC</th><th>BH P</th></tr></thead><tbody>"+
    pp.map(p=>{const r=statFor(p.protein_id);return `<tr><td>${esc(p.gene)}</td><td><button data-open-protein="${esc(p.protein_id)}" title="${esc(p.protein_id)}">${esc(p.protein_id)}</button></td><td>${fmt(r?.log2FC)}</td><td>${fmt(r?.["adj.P.Val"])}</td></tr>`;}).join("")+"</tbody></table>"+
    (!pp.length?'<p class="empty-note">没有可映射的唯一基因蛋白。</p>':"");
}
function displayTable(container,rows,fields,key=null,module=null){
  el(container).innerHTML="<table><thead><tr>"+fields.map(([k,label])=>"<th>"+esc(label)+"</th>").join("")+"</tr></thead><tbody>"+
    rows.slice(0,100).map(r=>"<tr>"+fields.map(([k])=>"<td>"+(k===key?
      `<button class="term-button" data-module="${module}" data-term="${esc(r[k])}">${esc(r[k])}</button>`:
      (finite(r[k])?fmt(r[k]):esc(r[k]??"—")))+"</td>").join("")+"</tr>").join("")+"</tbody></table>"+
    `<p class="chart-foot">共 ${rows.length} 条；表内最多展示100条，导出包含当前筛选的全部结果。</p>`;
}
function bind(id,handler){const node=el(id);if(!node._reportLinked){node.on("plotly_click",handler);node._reportLinked=true;}}
async function empty(id,message){
  const layout=baseLayout();layout.xaxis.visible=false;layout.yaxis.visible=false;
  layout.annotations=[{xref:"paper",yref:"paper",x:.5,y:.5,text:esc(message),showarrow:false,font:{size:13,color:"#617285"}}];
  await Plotly.react(el(id),[],layout,config);
}
async function openProtein(id){await C.selectProtein(id);await C.showView("protein");}
async function openSample(id){await C.selectSample(id);await C.showView("sample");}
function saveTsv(name,rows){
  const fields=[...new Set(rows.flatMap(r=>Object.keys(r)))].filter(k=>k!=="hits");
  const clean=v=>{if(v==null)return "NA";if(typeof v==="number")return String(v);let s=String(v).replace(/[\t\r\n]/g," ");return /^[=+\-@]/.test(s)?"'"+s:s;};
  const text=fields.join("\t")+"\n"+rows.map(r=>fields.map(k=>clean(r[k])).join("\t")).join("\n");
  const url=URL.createObjectURL(new Blob(["\uFEFF"+text],{type:"text/tab-separated-values;charset=utf-8"}));
  const a=document.createElement("a");a.href=url;a.download=name+".tsv";a.click();setTimeout(()=>URL.revokeObjectURL(url),2000);
}
async function qc(){
  const metric=el("qc-metric").value,layout=baseLayout(),ss=D.samples;
  layout.xaxis.tickangle=-40;layout.yaxis.title={text:metric==="observed"?"检出蛋白数":"缺失比例"};
  if(metric==="missing_fraction"){layout.yaxis.tickformat=".0%";layout.yaxis.range=[0,1];}
  await Plotly.react(el("qc-bars"),[{type:"bar",x:ss.map(s=>s.sample),y:ss.map(s=>s[metric]),
    customdata:ss.map(s=>s.sample),marker:{color:ss.map(s=>groupColor(s.group))},
    hovertemplate:"%{x}<br>%{y}<extra></extra>"}],layout,config);
  bind("qc-bars",e=>{if(e.points?.[0]?.customdata)openSample(e.points[0].customdata);});
  const traces=ss.map(s=>({type:"box",name:s.sample,y:Object.values(D.abundance).map(v=>v[s.sample]).filter(finite),
    boxpoints:false,marker:{color:groupColor(s.group)},line:{width:1.2},showlegend:false}));
  const bl=baseLayout();bl.xaxis.tickangle=-40;bl.yaxis.title={text:"已处理的 log2 丰度"};
  await Plotly.react(el("qc-boxes"),traces,bl,config);
  const names=ss.map(s=>s.sample),map=new Map(E.correlation.map(r=>[r[Object.keys(r)[0]],r]));
  const z=names.map(a=>names.map(b=>map.get(a)?.[b]??null));
  const cl=baseLayout();cl.xaxis.tickangle=-45;cl.yaxis.autorange="reversed";
  await Plotly.react(el("qc-correlation"),[{type:"heatmap",x:names,y:names,z,zmin:-1,zmax:1,colorscale:diverging,
    colorbar:{title:{text:"Pearson r"}},hoverongaps:false,hovertemplate:"%{y} × %{x}<br>r = %{z:.3f}<extra></extra>"}],cl,config);
  bind("qc-correlation",e=>{if(e.points?.[0]?.x)openSample(e.points[0].x);});
}
async function heatmap(){
  const mode=el("hm-set").value,q=el("hm-search").value.trim().toLowerCase().split(/[,，]/).filter(Boolean),count=+el("hm-count").value,scale=el("hm-scale").value;
  let rr=stats().filter(r=>mode==="sig"?["Up","Down"].includes(r.direction):variance(Object.values(D.abundance[r.protein_id]))>0);
  if(q.length)rr=rr.filter(r=>{const p=proteinMap.get(r.protein_id);return q.some(q=>`${p.protein_id} ${p.gene}`.toLowerCase().includes(q));});
  rr.sort((a,b)=>mode==="sig"?(a["adj.P.Val"]??1)-(b["adj.P.Val"]??1):variance(Object.values(D.abundance[b.protein_id]))-variance(Object.values(D.abundance[a.protein_id])));
  const total=rr.length;rr=rr.slice(0,count);
  el("hm-note").textContent=`${mode==="sig"?"正式差异蛋白，按BH校正P排序":"高变蛋白，按观测丰度方差排序；这是探索性展示"}。匹配 ${total} 个，展示 ${rr.length} 个。${scale==="z"?"Z-score采用观测值均值和样本标准差；颜色范围−2至2，超出范围颜色饱和。":"显示实际log2丰度。"}空白为缺失。`;
  if(!rr.length){await empty("protein-heatmap","当前条件下没有蛋白；可切换集合或调整关键词。");return;}
  const sampleIds=D.samples.map(s=>s.sample),z=[],custom=[];
  for(const r of rr){
    const v=sampleIds.map(s=>D.abundance[r.protein_id][s]),obs=v.filter(finite),mean=obs.reduce((s,x)=>s+x,0)/obs.length,sd=Math.sqrt(variance(obs));
    z.push(v.map(n=>!finite(n)?null:scale==="z"?(sd>0?(n-mean)/sd:0):n));
    custom.push(v.map(n=>[r.protein_id,n]));
  }
  const layout=baseLayout();layout.margin.l=185;layout.xaxis.tickangle=-45;layout.yaxis.autorange="reversed";
  fitLabels(layout,"protein-heatmap",rr.map(r=>(proteinMap.get(r.protein_id).gene||r.protein_id)+" · "+r.protein_id));
  layout.height=Math.max(420,rr.length*17+160);el("protein-heatmap").style.height=layout.height+"px";
  await Plotly.react(el("protein-heatmap"),[{type:"heatmap",x:sampleIds,y:rr.map(r=>(proteinMap.get(r.protein_id).gene||r.protein_id)+" · "+r.protein_id),
    z,customdata:custom,colorscale:scale==="z"?diverging:sequential,zmin:scale==="z"?-2:undefined,zmax:scale==="z"?2:undefined,
    hoverongaps:false,colorbar:{title:{text:scale==="z"?"行 Z":"log2"}},
    hovertemplate:"%{y}<br>样本: %{x}<br>显示值: %{z:.3f}<br>log2丰度: %{customdata[1]:.3f}<extra></extra>"}],layout,config);
  bind("protein-heatmap",e=>{const id=e.points?.[0]?.customdata?.[0];if(id)openProtein(id);});
}
function oraRows(){
  const db=el("ora-db").value,dir=el("ora-direction").value,q=el("ora-search").value.toLowerCase();
  return (E.ora[db]?.[coName()]?.[dir]||[]).filter(r=>r.overlap>0&&r.term.toLowerCase().includes(q)&&
    (el("ora-significance").value==="all"||(finite(r.padj)&&r.padj<=D.alpha))).slice().sort((a,b)=>(a.padj??1)-(b.padj??1));
}
function oraDetail(row){
  S.ora=row?.term??null;el("ora-title").textContent=row?.term||"当前没有符合条件的条目";
  el("ora-stats").innerHTML=row?smallStat("命中基因",row.overlap)+smallStat("背景基因",row.background_genes)+smallStat("BH P",fmt(row.padj)):"";
  el("ora-proteins").innerHTML=row?proteinTable(toGenes(row.genes)):"";
}
async function ora(){
  const rr=oraRows(),show=rr.slice(0,20);
  el("ora-note").textContent=`筛选后 ${rr.length} 条，图中展示BH校正P最小的前20条。横轴为命中数/本次选中且可注释的差异基因数；大小为命中数，颜色为−log10 BH校正P。`;
  displayTable("ora-table",rr,[["term","条目"],["overlap","命中数"],["term_background","条目背景"],["padj","BH P"]],"term","ora");
  oraDetail(rr.find(r=>r.term===S.ora)||rr[0]);
  if(!show.length){await empty("ora-plot","没有符合筛选条件的条目；无正式差异时不进行该富集检验。");return;}
  const layout=baseLayout();layout.margin.l=230;layout.yaxis.autorange="reversed";layout.xaxis.title={text:"Gene ratio"};layout.yaxis.tickfont={size:10};
  fitLabels(layout,"ora-plot",show.map(r=>r.term));
  await Plotly.react(el("ora-plot"),[{type:"scatter",mode:"markers",x:show.map(r=>r.overlap/r.selected_genes),y:show.map(r=>r.term),
    customdata:show.map(r=>r.term),text:show.map(r=>`${esc(r.term)}<br>命中 ${r.overlap}/${r.selected_genes}<br>BH P ${fmt(r.padj)}`),
    marker:{size:show.map(r=>10+Math.sqrt(r.overlap)*3),color:show.map(r=>-Math.log10(Math.max(r.padj,1e-300))),colorscale:sequential,showscale:true,colorbar:{title:{text:"−log10 BH P"}}},
    hovertemplate:"%{text}<extra></extra>"}],layout,config);
  bind("ora-plot",e=>oraDetail(oraRows().find(r=>r.term===e.points?.[0]?.customdata)));
}
const gseaSet=()=>E.gsea[el("gsea-db").value]?.[coName()];
function gseaRows(){
  const q=el("gsea-search").value.toLowerCase();
  return (gseaSet()?.rows||[]).filter(r=>r.pathway.toLowerCase().includes(q)&&(el("gsea-significance").value==="all"||(finite(r.padj)&&r.padj<=D.alpha)))
    .slice().sort((a,b)=>(a.padj??Infinity)-(b.padj??Infinity));
}
async function gseaDetail(row){
  S.gsea=row?.pathway??null;el("gsea-title").textContent=row?.pathway||"当前没有符合条件的通路";
  el("gsea-stats").innerHTML=row?smallStat("NES",fmt(row.NES))+smallStat("BH P",fmt(row.padj))+smallStat("集合大小",row.size):"";
  el("gsea-leading").innerHTML=row?proteinTable(toGenes(row.leadingEdge)):"";el("gsea-gene-focus").innerHTML="";
  if(!row){await empty("gsea-curve","选择通路后显示标准富集曲线");return;}
  const set=gseaSet(),n=set.values.length,hits=new Set(row.hits),total=row.hits.reduce((s,i)=>s+Math.abs(set.values[i]),0),nh=hits.size;
  let running=0;const curve=[0],rank=[0],hitx=[],hity=[],hitgenes=[];
  for(let i=0;i<n;i++){
    running+=hits.has(i)?(total?Math.abs(set.values[i])/total:1/nh):-1/(n-nh);
    curve.push(running);rank.push(i+1);
    if(hits.has(i)){hitx.push(i+1);hity.push(.5);hitgenes.push(set.genes[i]);}
  }
  const layout=baseLayout();layout.margin={l:70,r:25,t:25,b:55};layout.showlegend=false;
  layout.xaxis={...axis,domain:[0,1],anchor:"y3",title:{text:"基因排序位置（按有方向 t 值递减）"},range:[0,n+1]};
  layout.yaxis={...axis,domain:[.43,1],title:{text:"运行 ES"},zeroline:true};
  layout.yaxis2={visible:false,domain:[.30,.39],range:[0,1]};
  layout.yaxis3={...axis,domain:[0,.23],title:{text:"t 值"},zeroline:true};
  await Plotly.react(el("gsea-curve"),[
    {type:"scatter",mode:"lines",x:rank,y:curve,line:{color:"#009E73",width:2},hovertemplate:"排序 %{x}<br>ES %{y:.4f}<extra></extra>"},
    {type:"scatter",mode:"markers",x:hitx,y:hity,yaxis:"y2",customdata:hitgenes,marker:{symbol:"line-ns-open",size:12,color:"#34475a"},
      text:hitgenes.map(esc),hovertemplate:"%{text}<br>排序 %{x}<extra>点击查看蛋白</extra>"},
    {type:"scatter",mode:"lines",x:set.values.map((_,i)=>i+1),y:set.values,yaxis:"y3",fill:"tozeroy",line:{color:"#769cb4",width:1},
      hovertemplate:"排序 %{x}<br>t %{y:.3f}<extra></extra>"}],layout,config);
  bind("gsea-curve",e=>{const gene=e.points?.[0]?.customdata;if(gene)el("gsea-gene-focus").innerHTML="<h3>"+esc(gene)+"</h3>"+proteinTable([gene]);});
}
async function gsea(){
  const rr=gseaRows(),show=rr.filter(r=>finite(r.NES)).slice(0,20);
  el("gsea-note").textContent=`筛选后 ${rr.length} 条，展示前20条；原始NES与BH校正P直接来自fgsea。选择“全部受检条目”时包含未显著结果。`;
  displayTable("gsea-table",rr,[["pathway","通路"],["NES","NES"],["padj","BH P"],["size","基因数"]],"pathway","gsea");
  await gseaDetail(rr.find(r=>r.pathway===S.gsea)||rr[0]);
  if(!show.length){await empty("gsea-overview","没有符合筛选条件的通路，可切换显示范围。");return;}
  const layout=baseLayout();layout.margin.l=230;layout.yaxis.autorange="reversed";layout.xaxis.title={text:"NES"};layout.yaxis.tickfont={size:10};
  fitLabels(layout,"gsea-overview",show.map(r=>r.pathway));
  await Plotly.react(el("gsea-overview"),[{type:"bar",orientation:"h",x:show.map(r=>r.NES),y:show.map(r=>r.pathway),
    customdata:show.map(r=>r.pathway),marker:{color:show.map(r=>r.NES>=0?directionColor.Up:directionColor.Down)},
    text:show.map(r=>`${esc(r.pathway)}<br>NES ${fmt(r.NES)}<br>BH P ${fmt(r.padj)}`),hovertemplate:"%{text}<extra></extra>"}],layout,config);
  bind("gsea-overview",e=>gseaDetail(gseaRows().find(r=>r.pathway===e.points?.[0]?.customdata)));
}
const gsvaSet=()=>E.gsva[el("gsva-db").value];
function gsvaRows(){
  const set=gsvaSet();if(!set)return [];
  const diffs=new Map((set.differences[coName()]||[]).map(r=>[r.pathway,r]));
  const q=el("gsva-search").value.toLowerCase(),mode=el("gsva-order").value;
  return Object.keys(set.scores).map(pathway=>({...diffs.get(pathway),pathway,variance:variance(Object.values(set.scores[pathway]))}))
    .filter(r=>r.pathway.toLowerCase().includes(q)&&(mode!=="significant"||(finite(r["adj.P.Val"])&&r["adj.P.Val"]<=D.alpha)))
    .sort((a,b)=>mode==="variance"?b.variance-a.variance:(a["adj.P.Val"]??Infinity)-(b["adj.P.Val"]??Infinity));
}
async function gsvaDetail(row){
  S.gsva=row?.pathway??null;el("gsva-title").textContent=row?.pathway||"当前没有符合条件的通路";
  el("gsva-stats").innerHTML=row?smallStat("处理−对照分数差",fmt(row.score_difference))+smallStat("BH P",fmt(row["adj.P.Val"])):"";
  el("gsva-values").innerHTML="";
  if(!row){await empty("gsva-abundance","选择通路后显示各样本分数");return;}
  const values=gsvaSet().scores[row.pathway],co=contrast(),ss=D.samples.filter(s=>[co.control,co.treatment].includes(s.group)),traces=[];
  if(co.paired)for(const subject of new Set(ss.map(s=>s.subject))){
    const a=ss.find(s=>s.subject===subject&&s.group===co.control),b=ss.find(s=>s.subject===subject&&s.group===co.treatment);
    if(a&&b&&finite(values[a.sample])&&finite(values[b.sample]))traces.push({type:"scatter",mode:"lines",x:[co.control,co.treatment],y:[values[a.sample],values[b.sample]],line:{color:"#c6d0da",width:1.4},hoverinfo:"skip",showlegend:false});
  }
  for(const g of [co.control,co.treatment]){
    const a=ss.filter(s=>s.group===g);
    traces.push({type:"scatter",mode:"markers",x:a.map(()=>g),y:a.map(s=>values[s.sample]),
      customdata:a.map(s=>s.sample),text:a.map(s=>esc(s.sample)),marker:{size:10,color:groupColor(g)},name:g,
      hovertemplate:"%{text}<br>GSVA %{y:.4f}<extra></extra>"});
  }
  const layout=baseLayout();layout.yaxis.title={text:"GSVA 分数"};layout.xaxis.categoryorder="array";layout.xaxis.categoryarray=[co.control,co.treatment];
  await Plotly.react(el("gsva-abundance"),traces,layout,config);
  bind("gsva-abundance",e=>{const id=e.points?.[0]?.customdata;if(id)openSample(id);});
  el("gsva-values").innerHTML="<table><thead><tr><th>样本</th><th>组别</th><th>GSVA分数</th></tr></thead><tbody>"+
    D.samples.map(s=>`<tr><td>${esc(s.sample)}</td><td>${esc(s.group)}</td><td>${fmt(values[s.sample])}</td></tr>`).join("")+"</tbody></table>";
}
async function gsva(){
  const rr=gsvaRows(),show=rr.slice(0,30),set=gsvaSet(),ss=D.samples.map(s=>s.sample);
  el("gsva-coverage").textContent=`完整观测筛选后 ${E.gsva_gene_count} 个基因；当前数据库共 ${Object.keys(set?.scores||{}).length} 条可评分通路。`;
  el("gsva-note").textContent=`筛选后 ${rr.length} 条，热图展示前30条。颜色为原始GSVA分数，不进行行Z-score转换；以0为色阶中心。按方差展示不等于显著差异。`;
  displayTable("gsva-table",rr,[["pathway","通路"],["score_difference","处理−对照分数差"],["adj.P.Val","BH P"],["variance","方差"]],"pathway","gsva");
  await gsvaDetail(rr.find(r=>r.pathway===S.gsva)||rr[0]);
  if(!show.length){await empty("gsva-heatmap","没有符合条件的通路，或此模块未配置。");return;}
  const z=show.map(r=>ss.map(s=>set.scores[r.pathway][s])),max=Math.max(.01,...z.flat().filter(finite).map(Math.abs));
  const layout=baseLayout();layout.margin.l=260;layout.xaxis.tickangle=-45;layout.yaxis.autorange="reversed";layout.yaxis.tickfont={size:10};
  fitLabels(layout,"gsva-heatmap",show.map(r=>r.pathway));
  layout.height=Math.max(430,show.length*19+150);el("gsva-heatmap").style.height=layout.height+"px";
  await Plotly.react(el("gsva-heatmap"),[{type:"heatmap",x:ss,y:show.map(r=>r.pathway),z,colorscale:diverging,zmin:-max,zmax:max,
    colorbar:{title:{text:"GSVA"}},hovertemplate:"%{y}<br>%{x}<br>GSVA %{z:.4f}<extra></extra>"}],layout,config);
  bind("gsva-heatmap",e=>gsvaDetail(gsvaRows().find(r=>r.pathway===e.points?.[0]?.y)));
}
async function ppi(){
  const set=E.ppi[coName()]||{nodes:[],edges:[]},nodes=set.nodes,edges=set.edges.slice().sort((a,b)=>b.score-a.score).slice(0,100);
  el("ppi-note").textContent=`完整结果 ${set.edges.length} 条边；图中显示最高分的至多100条边。线仅表示数据库中的关联，不证明本实验中的因果关系。`;
  displayTable("ppi-table",set.edges,[["gene1","基因1"],["gene2","基因2"],["score","STRING分数"]]);
  el("ppi-gene").textContent="点击节点查看对应蛋白";el("ppi-proteins").innerHTML="";
  if(!nodes.length){await empty("ppi-plot","没有满足阈值的正式差异基因相互作用，不生成网络。");return;}
  const nmap=new Map(nodes.map(n=>[n.gene,n])),xx=[],yy=[];
  for(const e of edges){const a=nmap.get(e.gene1),b=nmap.get(e.gene2);xx.push(a.x,b.x,null);yy.push(a.y,b.y,null);}
  const colors=nodes.map(n=>{
    const dir=new Set(geneProteins([n.gene]).map(p=>statFor(p.protein_id)?.direction).filter(d=>["Up","Down"].includes(d)));
    return dir.size===1?directionColor[[...dir][0]]:"#7b6ea8";
  });
  const layout=baseLayout();layout.xaxis.visible=false;layout.yaxis.visible=false;layout.showlegend=false;layout.margin={l:40,r:40,t:35,b:35};
  await Plotly.react(el("ppi-plot"),[
    {type:"scatter",mode:"lines",x:xx,y:yy,line:{color:"#c3ced8",width:1.2},hoverinfo:"skip"},
    {type:"scatter",mode:"markers+text",x:nodes.map(n=>n.x),y:nodes.map(n=>n.y),text:nodes.map(n=>esc(n.gene)),
      customdata:nodes.map(n=>n.gene),textposition:"top center",marker:{size:nodes.map(n=>13+Math.sqrt(n.degree)*3),color:colors,line:{color:"#fff",width:1.5}},
      hovertemplate:"%{text}<extra>点击查看蛋白</extra>"}],layout,config);
  bind("ppi-plot",e=>{const gene=e.points?.[0]?.customdata;if(gene){S.ppi=gene;el("ppi-gene").textContent=gene;el("ppi-proteins").innerHTML=proteinTable([gene]);}});
}
function files(){
  displayTable("summary-table",E.summary,[["comparison","比较"],["tested","受检"],["up","升高"],["down","降低"],["FDR_DEP","正式差异"],["raw_p_candidates","未校正候选"]]);
  const rows=[...E.statuses.function,...E.statuses.pathways];
  displayTable("module-status",rows,[["module","模块"],["database","数据库"],["comparison","比较"],["status","状态"],["tested_terms","受检条目"]]);
  const q=el("file-search").value.toLowerCase();
  el("file-list").innerHTML=[...E.files,...E.image_files].filter(f=>f.label.toLowerCase().includes(q))
    .map(f=>`<a href="${esc(f.path)}">${esc(f.label)}</a>`).join("");
}
const views={qc,heatmap,ora,gsea,gsva,ppi,files};
async function render(view){if(views[view])await views[view]();}
window.proteomicsExtensions={render,oraRows,gseaRows,gsvaRows,getSelection:()=>({...S})};
C.ready.then(()=>{
  for(const name of Object.keys(views))el("nav-"+name).addEventListener("click",()=>C.showView(name));
  document.querySelectorAll(".shared-contrast").forEach(s=>{
    s.innerHTML=D.contrasts.map(c=>`<option value="${esc(c.name)}">${esc(c.treatment)} / ${esc(c.control)}</option>`).join("");
    s.addEventListener("change",()=>C.setContrast(s.value));
  });
  document.querySelectorAll(".database-select").forEach(s=>{
    s.innerHTML=E.databases.map(db=>`<option value="${esc(db)}">${esc(db)}</option>`).join("");
  });
  const controls={qc:["qc-metric"],heatmap:["hm-set","hm-count","hm-scale","hm-search"],
    ora:["ora-db","ora-direction","ora-significance","ora-search"],gsea:["gsea-db","gsea-significance","gsea-search"],
    gsva:["gsva-db","gsva-order","gsva-search"],files:["file-search"]};
  for(const [view,ids] of Object.entries(controls))for(const id of ids)el(id).addEventListener(el(id).tagName==="INPUT"?"input":"change",()=>render(view));
  document.addEventListener("click",e=>{
    const protein=e.target.closest("[data-open-protein]");if(protein){openProtein(protein.dataset.openProtein);return;}
    const term=e.target.closest("[data-module]");if(!term)return;
    if(term.dataset.module==="ora")oraDetail(oraRows().find(r=>r.term===term.dataset.term));
    if(term.dataset.module==="gsea")gseaDetail(gseaRows().find(r=>r.pathway===term.dataset.term));
    if(term.dataset.module==="gsva")gsvaDetail(gsvaRows().find(r=>r.pathway===term.dataset.term));
  });
  el("ora-export").addEventListener("click",()=>saveTsv(coName()+"."+el("ora-db").value+".ORA",oraRows()));
  el("gsea-export").addEventListener("click",()=>saveTsv(coName()+"."+el("gsea-db").value+".GSEA",gseaRows()));
  el("gsva-export").addEventListener("click",()=>saveTsv(coName()+"."+el("gsva-db").value+".GSVA",gsvaRows()));
  el("ppi-export").addEventListener("click",()=>saveTsv(coName()+".PPI",E.ppi[coName()]?.edges||[]));
  el("gsea-curve-export").addEventListener("click",()=>Plotly.downloadImage(el("gsea-curve"),{format:"svg",filename:"GSEA_enrichment",width:1100,height:700}));
});
})();
