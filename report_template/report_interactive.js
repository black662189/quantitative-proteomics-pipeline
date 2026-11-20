/* Plotly-based view layer. Statistical results and PCA coordinates are read-only. */
(() => {
"use strict";
const data=JSON.parse(document.getElementById("report-data").textContent);
const el=id=>document.getElementById(id);
const esc=s=>String(s??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const finite=n=>typeof n==="number" && Number.isFinite(n);
const fmt=(n,d=3)=>finite(n)?(n!==0 && Math.abs(n)<.001?n.toExponential(2):n.toFixed(d)):"—";
const state={contrast:data.contrasts[0].name,protein:null,sample:data.samples[0].sample,page:0,view:"protein"};
const proteinMap=new Map(data.proteins.map(p=>[p.protein_id,p]));
const sampleMap=new Map(data.samples.map(s=>[s.sample,s]));
const groups=[...new Set(data.samples.map(s=>s.group))].sort();
const palette=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9"];
const groupColor=g=>palette[Math.max(0,groups.indexOf(g))%palette.length];
const directionColor={Up:"#D55E00",Down:"#0072B2",NS:"#B8C2CB",Not_tested:"#B8C2CB"};
const directionName={Up:"升高",Down:"降低",NS:"未显著",Not_tested:"未检验"};
const config={responsive:true,displaylogo:false,scrollZoom:false,modeBarButtonsToRemove:["lasso2d","select2d"],
  toImageButtonOptions:{format:"png",scale:2,filename:"proteomics_plot"}};
const axis={showline:true,linecolor:"#cbd5df",zeroline:false,gridcolor:"#edf1f5",ticks:"outside",tickcolor:"#cbd5df",automargin:true};
const baseLayout=()=>({font:{family:"Segoe UI, Microsoft YaHei, sans-serif",size:12,color:"#34475a"},
  paper_bgcolor:"#fff",plot_bgcolor:"#fff",margin:{l:68,r:24,t:22,b:62},
  hovermode:"closest",hoverlabel:{bgcolor:"#fff",font:{size:12}},
  legend:{orientation:"h",x:0,y:1.12,font:{size:11}},
  xaxis:{...axis},yaxis:{...axis}});
const contrast=()=>data.contrasts.find(c=>c.name===state.contrast);
const stats=()=>data.statistics[state.contrast];
const statFor=id=>stats().find(r=>r.protein_id===id);
function metric(label,value,color=""){return `<div class="metric"><strong${color?` style="color:${color}"`:""}>${esc(value)}</strong><span>${esc(label)}</span></div>`;}
function smallStat(label,value){return `<div><span>${esc(label)}</span><strong>${esc(value)}</strong></div>`;}
function pointValue(r){return finite(r.log2FC)&&finite(r["adj.P.Val"])&&r.tested;}
async function renderVolcano(){
  const co=contrast(),rows=stats(),traces=[];
  for(const dir of ["NS","Down","Up"]){
    const a=rows.filter(r=>pointValue(r)&&r.direction===dir);
    traces.push({type:"scatter",mode:"markers",name:directionName[dir],
      x:a.map(r=>r.log2FC),y:a.map(r=>-Math.log10(Math.max(r["adj.P.Val"],1e-300))),
      customdata:a.map(r=>r.protein_id),text:a.map(r=>{
        const p=proteinMap.get(r.protein_id);
        return `${esc(p.gene||p.protein_id)}<br>${esc(p.protein_id)}<br>log2FC: ${fmt(r.log2FC)}<br>BH P: ${fmt(r["adj.P.Val"])}`;
      }),hovertemplate:"%{text}<extra></extra>",marker:{size:5.5,color:directionColor[dir],opacity:dir==="NS"?.5:.9}});
  }
  const sel=statFor(state.protein);
  if(sel&&pointValue(sel))traces.push({type:"scatter",mode:"markers",name:"所选蛋白",showlegend:false,
    x:[sel.log2FC],y:[-Math.log10(Math.max(sel["adj.P.Val"],1e-300))],customdata:[sel.protein_id],
    marker:{size:13,color:"rgba(0,0,0,0)",line:{color:"#25374b",width:2}},
    text:[esc(sel.protein_id)],hovertemplate:"%{text}<extra>所选蛋白</extra>"});
  const cut=Math.log2(data.fc),ycut=-Math.log10(data.alpha);
  const layout=baseLayout();layout.xaxis.title={text:`log2FC（${esc(co.treatment)} / ${esc(co.control)}）`};
  layout.yaxis.title={text:"−log10（BH 校正 P）"};layout.uirevision=state.contrast;
  layout.shapes=[-cut,cut].map(v=>({type:"line",xref:"x",yref:"paper",x0:v,x1:v,y0:0,y1:1,line:{color:"#93a2b3",width:1,dash:"dot"}}));
  layout.shapes.push({type:"line",xref:"paper",yref:"y",x0:0,x1:1,y0:ycut,y1:ycut,line:{color:"#93a2b3",width:1,dash:"dot"}});
  el("comparison-note").textContent=`${co.treatment} 相对 ${co.control} · ${co.paired?"配对设计":"独立组设计"} · BH P ≤ ${data.alpha} 且 |log2FC| ≥ ${cut.toFixed(3)}`;
  const n=rows.filter(pointValue).length,up=rows.filter(r=>r.direction==="Up").length,down=rows.filter(r=>r.direction==="Down").length;
  el("contrast-metrics").innerHTML=metric("受检蛋白",n)+metric("升高",up,directionColor.Up)+metric("降低",down,directionColor.Down)+metric("未检验",rows.length-n);
  el("volcano-note").textContent=`绘制 ${n} 个受检蛋白；${rows.length-n} 个未检验蛋白保留在结果表中。搜索和筛选仅作用于下表，不重新计算统计结果。`;
  await Plotly.react(el("volcano"),traces,layout,config);
}
async function renderProtein(){
  const p=proteinMap.get(state.protein),r=statFor(state.protein),co=contrast();
  if(!p||!r)return;
  el("protein-gene").textContent=p.gene||"无唯一基因注释";
  el("protein-id").textContent=p.protein_id;el("protein-description").textContent=p.description||"未提供蛋白描述";
  el("protein-direction").textContent=directionName[r.direction]||r.direction;
  el("protein-direction").className="badge "+(r.direction==="Up"?"up":r.direction==="Down"?"down":"");
  el("protein-statistics").innerHTML=smallStat("log2FC",fmt(r.log2FC))+smallStat("原始 P",fmt(r["P.Value"]))+
    smallStat("BH 校正 P",fmt(r["adj.P.Val"]))+(co.paired?smallStat("完整配对",r.n_complete_pairs??"—"):"");
  const ss=data.samples.filter(s=>[co.control,co.treatment].includes(s.group));
  const values=data.abundance[p.protein_id],pos=new Map(),traces=[];
  for(const [i,g] of [co.control,co.treatment].entries()){
    const a=ss.filter(s=>s.group===g).sort((a,b)=>a.sample.localeCompare(b.sample));
    a.forEach((s,k)=>pos.set(s.sample,i+(a.length===1?0:(k/(a.length-1)-.5)*.16)));
  }
  if(co.paired){
    for(const subject of new Set(ss.map(s=>s.subject))){
      const a=ss.find(s=>s.subject===subject&&s.group===co.control),b=ss.find(s=>s.subject===subject&&s.group===co.treatment);
      if(a&&b&&finite(values[a.sample])&&finite(values[b.sample]))traces.push({type:"scatter",mode:"lines",showlegend:false,hoverinfo:"skip",
        x:[pos.get(a.sample),pos.get(b.sample)],y:[values[a.sample],values[b.sample]],line:{color:"#c6d0da",width:1.4}});
    }
  }
  for(const g of [co.control,co.treatment]){
    const a=ss.filter(s=>s.group===g&&finite(values[s.sample]));
    traces.push({type:"scatter",mode:"markers",name:g,x:a.map(s=>pos.get(s.sample)),y:a.map(s=>values[s.sample]),
      customdata:a.map(s=>s.sample),text:a.map(s=>`${esc(s.sample)}<br>组别: ${esc(g)}<br>个体: ${esc(s.subject||"—")}<br>log2丰度: ${fmt(values[s.sample])}`),
      hovertemplate:"%{text}<extra></extra>",marker:{color:groupColor(g),size:10,line:{color:"#fff",width:1.3}}});
  }
  const layout=baseLayout();layout.margin={l:65,r:18,t:35,b:42};
  layout.xaxis={...axis,tickmode:"array",tickvals:[0,1],ticktext:[co.control,co.treatment],range:[-.4,1.4],showgrid:false};
  layout.yaxis.title={text:"已处理的 log2 丰度"};layout.showlegend=false;
  if(!ss.some(s=>finite(values[s.sample])))layout.annotations=[{text:"当前两组均无观测丰度",xref:"paper",yref:"paper",x:.5,y:.5,showarrow:false}];
  await Plotly.react(el("protein-abundance"),traces,layout,config);
  el("protein-note").textContent=`仅展示当前比较的两组，${ss.filter(s=>!finite(values[s.sample])).length} 个样本缺失；不填补。${co.paired?"灰线连接同一个体。":""}点击丰度点可转到样本详情。`;
  el("abundance-table").innerHTML="<table><thead><tr><th>样本</th><th>组别</th><th>log2 丰度</th></tr></thead><tbody>"+
    data.samples.map(s=>`<tr><td>${esc(s.sample)}</td><td>${esc(s.group)}</td><td>${finite(values[s.sample])?fmt(values[s.sample]):"缺失"}</td></tr>`).join("")+"</tbody></table>";
}
function filteredRows(){
  const q=el("protein-search").value.trim().toLowerCase(),filter=el("protein-filter").value,sort=el("protein-sort").value;
  return stats().filter(r=>{
    const p=proteinMap.get(r.protein_id);
    return (!q||`${p.protein_id} ${p.gene} ${p.gene_annotation} ${p.description}`.toLowerCase().includes(q))&&
      (filter==="all"||(filter==="sig"?["Up","Down"].includes(r.direction):r.direction===filter));
  }).slice().sort((a,b)=>{
    if(sort==="gene")return (proteinMap.get(a.protein_id).gene||a.protein_id).localeCompare(proteinMap.get(b.protein_id).gene||b.protein_id);
    if(sort==="effect")return (finite(b.log2FC)?Math.abs(b.log2FC):-1)-(finite(a.log2FC)?Math.abs(a.log2FC):-1);
    return (a["adj.P.Val"]??Infinity)-(b["adj.P.Val"]??Infinity)||a.protein_id.localeCompare(b.protein_id);
  });
}
function renderTable(){
  const rows=filteredRows(),pageSize=15,pages=Math.max(1,Math.ceil(rows.length/pageSize));
  state.page=Math.min(state.page,pages-1);
  el("table-status").textContent=`匹配 ${rows.length} / ${stats().length} 个蛋白。点击蛋白 ID 查看详情。`;
  el("protein-rows").innerHTML=rows.slice(state.page*pageSize,(state.page+1)*pageSize).map(r=>
    `<tr class="${r.protein_id===state.protein?"selected":""}"><td><button data-protein="${esc(r.protein_id)}" title="${esc(r.protein_id)}">${esc(r.protein_id)}</button></td><td>${esc(proteinMap.get(r.protein_id).gene||"—")}</td><td>${fmt(r.log2FC)}</td><td>${fmt(r["P.Value"])}</td><td>${fmt(r["adj.P.Val"])}</td><td>${esc(directionName[r.direction])}</td></tr>`).join("")||
    '<tr><td colspan="6">没有匹配结果，请调整关键词或筛选条件。</td></tr>';
  el("page-number").textContent=`${state.page+1} / ${pages}`;el("previous-page").disabled=state.page===0;el("next-page").disabled=state.page>=pages-1;
}
async function selectProtein(id){
  if(!proteinMap.has(id))return;
  state.protein=id;renderTable();await Promise.all([renderVolcano(),renderProtein()]);
}
async function renderPca(){
  const traces=[];
  if(!data.pca.length){
    const layout=baseLayout();layout.annotations=[{text:"完整观测的非恒定蛋白不足，未计算PCA",xref:"paper",yref:"paper",x:.5,y:.5,showarrow:false}];
    await Plotly.react(el("pca"),[],layout,config);el("pca-note").textContent="仍可从下拉框或样本表查看样本详情。";return;
  }
  for(const g of groups){
    const rows=data.pca.filter(p=>sampleMap.get(p.sample).group===g);
    traces.push({type:"scatter",mode:"markers",name:g,x:rows.map(p=>p.PC1),y:rows.map(p=>p.PC2),
      customdata:rows.map(p=>p.sample),text:rows.map(p=>`${esc(p.sample)}<br>组别: ${esc(g)}<br>PC1: ${fmt(p.PC1)}<br>PC2: ${fmt(p.PC2)}`),
      hovertemplate:"%{text}<extra></extra>",marker:{color:groupColor(g),size:12,line:{color:"#fff",width:1.5}}});
  }
  const s=data.pca.find(p=>p.sample===state.sample);
  traces.push({type:"scatter",mode:"markers",showlegend:false,x:[s.PC1],y:[s.PC2],customdata:[s.sample],
    marker:{symbol:"diamond-open",color:"#25374b",size:20,line:{width:2}},text:[esc(s.sample)],hovertemplate:"%{text}<extra>当前样本</extra>"});
  const layout=baseLayout();layout.xaxis.title={text:`PC1 (${(100*data.explained[0]).toFixed(1)}%)`};layout.yaxis.title={text:`PC2 (${(100*data.explained[1]).toFixed(1)}%)`};
  layout.uirevision="pca";await Plotly.react(el("pca"),traces,layout,config);
  el("pca-note").textContent=`使用 ${data.complete_proteins} 个完整观测且非恒定蛋白；按蛋白中心化、不做方差标准化、不填补。显示既有 PCA 坐标，不以组间分离程度判断显著性。`;
}
async function renderSample(){
  const s=sampleMap.get(state.sample);el("sample-select").value=s.sample;
  el("sample-name").textContent=s.sample;el("sample-group").textContent="组别 "+s.group;
  el("sample-metadata").innerHTML=[["组别",s.group],["个体 ID",s.subject||"未提供"],["批次 / plex",s.plex],["TMT 通道",s.channel||"不适用"]].map(([k,v])=>smallStat(k,v)).join("");
  el("sample-metrics").innerHTML=smallStat("检出蛋白",s.observed)+smallStat("缺失率",(s.missing_fraction*100).toFixed(1)+"%")+smallStat("log2 中位数",fmt(s.median_log2,2));
  const values=Object.values(data.abundance).map(row=>row[s.sample]).filter(finite);
  const layout=baseLayout();layout.margin={l:60,r:18,t:15,b:48};layout.xaxis.title={text:"已处理的 log2 丰度"};layout.yaxis.title={text:"蛋白数"};layout.showlegend=false;
  await Plotly.react(el("sample-distribution"),[{type:"histogram",x:values,nbinsx:30,
    marker:{color:groupColor(s.group),opacity:.8},hovertemplate:"log2 丰度区间: %{x}<br>蛋白数: %{y}<extra></extra>"}],layout,config);
  el("sample-table").innerHTML="<table><thead><tr><th>样本</th><th>组别</th><th>个体 ID</th><th>检出蛋白</th><th>缺失率</th><th>log2 中位数</th></tr></thead><tbody>"+
    data.samples.map(a=>`<tr class="${a.sample===state.sample?"selected":""}"><td><button data-sample="${esc(a.sample)}">${esc(a.sample)}</button></td><td>${esc(a.group)}</td><td>${esc(a.subject||"—")}</td><td>${a.observed}</td><td>${(100*a.missing_fraction).toFixed(1)}%</td><td>${fmt(a.median_log2,2)}</td></tr>`).join("")+"</tbody></table>";
}
async function selectSample(id){
  if(!sampleMap.has(id))return;
  state.sample=id;await Promise.all([renderPca(),renderSample()]);
}
async function showView(view){
  state.view=view;
  for(const name of ["protein","sample","qc","heatmap","ora","gsea","gsva","ppi","files"]){
    el(name+"-page").hidden=name!==view;el("nav-"+name).classList.toggle("active",name===view);el("nav-"+name).setAttribute("aria-pressed",String(name===view));
  }
  await new Promise(resolve=>requestAnimationFrame(resolve));
  if(!["sample","protein"].includes(view)){await window.proteomicsExtensions.render(view);return;}
  if(view==="sample")await selectSample(state.sample);else await selectProtein(state.protein);
  for(const id of view==="sample"?["pca","sample-distribution"]:["volcano","protein-abundance"])Plotly.Plots.resize(el(id));
}
async function setContrast(name){
  state.contrast=name;state.page=0;el("contrast-select").value=name;
  document.querySelectorAll(".shared-contrast").forEach(e=>e.value=name);
  await selectProtein(state.protein);
  if(window.proteomicsExtensions&&!["protein","sample"].includes(state.view))await window.proteomicsExtensions.render(state.view);
}
function downloadResults(){
  const rows=filteredRows(),fields=["protein_id","gene","log2FC","P.Value","adj.P.Val","direction","tested"];
  const clean=v=>{if(v==null)return "NA";if(typeof v==="number")return String(v);let s=String(v).replace(/[\t\r\n]/g," ");return /^[=+\-@]/.test(s)?"'"+s:s;};
  const text=fields.join("\t")+"\n"+rows.map(r=>fields.map(k=>clean(k==="gene"?proteinMap.get(r.protein_id).gene:r[k])).join("\t")).join("\n");
  const url=URL.createObjectURL(new Blob(["\uFEFF"+text],{type:"text/tab-separated-values;charset=utf-8"}));
  const a=document.createElement("a");a.href=url;a.download=state.contrast+".filtered.tsv";a.click();setTimeout(()=>URL.revokeObjectURL(url),2000);
}
async function boot(){
  el("project-title").textContent=data.project;
  const names={hsa:"人",fca:"猫",mmu:"小鼠",rno:"大鼠",cfa:"犬",ssc:"猪",mcf:"食蟹猴",mcc:"恒河猴"};
  el("species").textContent=`${names[data.species]||data.species} · TaxID ${data.taxid}`;
  el("project-counts").textContent=`${data.samples.length} 个样本 · ${data.proteins.length} 个保留蛋白`;
  el("contrast-select").innerHTML=data.contrasts.map(c=>`<option value="${esc(c.name)}">${esc(c.treatment)} / ${esc(c.control)}</option>`).join("");
  el("sample-select").innerHTML=data.samples.map(s=>`<option value="${esc(s.sample)}">${esc(s.sample)}</option>`).join("");
  el("processing-note").textContent=`输入尺度：${data.input_scale}；归一化设置：${data.normalization}。蛋白丰度图显示预处理后的 log2 值，不是原始线性强度。`;
  el("input-warnings").innerHTML=data.warnings.map(s=>"<li>"+esc(s)+"</li>").join("");
  state.protein=(stats().find(r=>r.direction==="Up"||r.direction==="Down")||stats()[0]).protein_id;
  await selectProtein(state.protein);
  el("volcano").on("plotly_click",e=>{const id=e.points?.[0]?.customdata;if(id)selectProtein(id);});
  el("protein-abundance").on("plotly_click",e=>{const id=e.points?.[0]?.customdata;if(id){state.sample=id;showView("sample");}});
  await selectSample(state.sample);
  el("pca").on("plotly_click",e=>{const id=e.points?.[0]?.customdata;if(id)selectSample(id);});
  el("nav-protein").addEventListener("click",()=>showView("protein"));
  el("nav-sample").addEventListener("click",()=>showView("sample"));
  el("contrast-select").addEventListener("change",e=>setContrast(e.target.value));
  el("sample-select").addEventListener("change",e=>selectSample(e.target.value));
  for(const id of ["protein-search","protein-filter","protein-sort"])el(id).addEventListener(id==="protein-search"?"input":"change",()=>{state.page=0;renderTable();});
  el("protein-rows").addEventListener("click",e=>{const b=e.target.closest("[data-protein]");if(b)selectProtein(b.dataset.protein);});
  el("sample-table").addEventListener("click",e=>{const b=e.target.closest("[data-sample]");if(b)selectSample(b.dataset.sample);});
  el("previous-page").addEventListener("click",()=>{state.page=Math.max(0,state.page-1);renderTable();});
  el("next-page").addEventListener("click",()=>{state.page++;renderTable();});
  el("download-results").addEventListener("click",downloadResults);
  el("export-volcano").addEventListener("click",()=>Plotly.downloadImage(el("volcano"),{format:"svg",filename:state.contrast+".volcano",width:1000,height:700}));
  el("export-pca").addEventListener("click",()=>Plotly.downloadImage(el("pca"),{format:"svg",filename:"PCA",width:1000,height:700}));
}
window.proteomicsPreview={ready:boot(),getState:()=>({...state}),filteredRows,
  data,el,esc,finite,fmt,contrast,stats,statFor,proteinMap,sampleMap,groupColor,directionColor,config,axis,baseLayout,smallStat,
  selectProtein,selectSample,showView,setContrast};
window.proteomicsPreview.ready.catch(err=>{console.error(err);el("comparison-note").textContent="页面初始化失败："+err.message;el("comparison-note").setAttribute("role","alert");});
})();
