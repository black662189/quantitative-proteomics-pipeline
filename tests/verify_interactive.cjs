// DOM and event-linkage tests; this does not claim real-browser rendering QA.
const fs=require("node:fs"),assert=require("node:assert/strict");
const {JSDOM}=require("jsdom");
const path=require("node:path");
const root=process.argv[2];
if(!root)throw Error("Usage: node tests/verify_interactive.cjs <directory containing TMT_HepG2 and DIA_cat>");
const pipeline=path.resolve(__dirname,"..");
const code=fs.readFileSync(path.join(pipeline,"report_template/report_interactive.js"),"utf8");
const modules=fs.readFileSync(path.join(pipeline,"report_template/report_modules.js"),"utf8");
// Load the real library only for schema validation; no rendering surface is used.
const schemaWindow=new JSDOM("",{runScripts:"outside-only"}).window;
schemaWindow.HTMLCanvasElement.prototype.getContext=()=>({measureText:s=>({width:String(s).length*7})});
schemaWindow.URL.createObjectURL=()=>"";
schemaWindow.eval(fs.readFileSync(path.join(pipeline,"vendor/plotly-3.1.0.min.js"),"utf8"));
const tick=()=>new Promise(resolve=>setTimeout(resolve,50));
const checkedSchemas=new Set();
async function check(name){
  const html=fs.readFileSync(`${root}/${name}/proteomics_report.html`,"utf8");
  assert(html.includes(code)&&html.includes(modules),"Embedded report scripts match deployed source");
  const dom=new JSDOM(html,{runScripts:"outside-only",pretendToBeVisual:true});
  const w=dom.window,d=w.document,calls=[],listeners=new Map();
  // Plotly is mocked at its documented API boundary so callbacks can be checked
  // without accessing a browser surface blocked from local files.
  w.Plotly={
    react:async(node,traces,layout,config)=>{
      const t=schemaWindow.JSON.parse(JSON.stringify(traces)),l=schemaWindow.JSON.parse(JSON.stringify(layout));
      const nonempty=t.filter(r=>r.x?.length||r.y?.length||r.z?.length);
      const schemaKey=node.id+":"+nonempty.map(r=>r.type+"/"+r.mode+"/"+r.yaxis).join(",")+":"+Object.keys(l).join(",");
      if(nonempty.length&&!checkedSchemas.has(schemaKey)){
        // Plotly.validate 3.1.0 crashes on nested attributes of suppressed
        // legends (also auto-suppressed single traces). Omit display-only
        // legend styling and hidden axes from this partial schema check.
        delete l.legend;
        for(const k of Object.keys(l))if(/^xaxis\d*$|^yaxis\d*$/.test(k)&&l[k].visible===false)l[k]={visible:false};
        let problems;
        try {problems=schemaWindow.Plotly.validate(nonempty,l);}
        catch(err){fs.writeFileSync(`${root}/schema_error.json`,JSON.stringify({plot:node.id,traces:nonempty,layout:l}));throw Error("Plotly schema validation failed for "+node.id+": "+err.message);}
        const errors=(problems||[]).filter(p=>p.code!=="unused");
        assert(!errors.length,JSON.stringify({plot:node.id,errors}));
        checkedSchemas.add(schemaKey);
      }
      node._traces=traces;node._layout=layout;node._config=config;
      node.on=(event,fn)=>listeners.set(node.id+":"+event,fn);
    },
    Plots:{resize:()=>{}},
    downloadImage:(node,options)=>{calls.push({node:node.id,options});return Promise.resolve();}
  };
  w.eval(code);w.eval(modules);await w.proteomicsPreview.ready;await tick();
  const payload=JSON.parse(d.getElementById("report-data").textContent);
  const get=id=>d.getElementById(id);
  const change=(id,value,type="change")=>{get(id).value=value;get(id).dispatchEvent(new w.Event(type,{bubbles:true}));};
  const click=id=>get(id).click();
  const current=()=>w.proteomicsPreview.getState();
  assert.equal(get("project-title").textContent,payload.project);
  const co=payload.contrasts[0],rows=payload.statistics[co.name];
  const tested=rows.filter(r=>r.tested&&r.log2FC!==null&&r["adj.P.Val"]!==null);
  assert.equal(get("volcano")._traces.slice(0,3).reduce((n,t)=>n+t.x.length,0),tested.length);
  const target=tested.at(-1);
  listeners.get("volcano:plotly_click")({points:[{customdata:target.protein_id}]});await tick();
  assert.equal(current().protein,target.protein_id);
  assert.equal(get("protein-id").textContent,target.protein_id);
  const inputAbundance=payload.abundance[target.protein_id];
  for(const trace of get("protein-abundance")._traces.filter(t=>t.mode==="markers")){
    trace.customdata.forEach((sample,i)=>assert.equal(trace.y[i],inputAbundance[sample]));
  }
  change("protein-search",target.protein_id,"input");assert(w.proteomicsPreview.filteredRows().some(r=>r.protein_id===target.protein_id));
  change("protein-search","unlikely_no_protein_123xyz","input");assert.equal(w.proteomicsPreview.filteredRows().length,0);
  assert(get("protein-rows").textContent.includes("没有匹配"));
  change("protein-search","","input");
  change("protein-filter","sig");assert.equal(w.proteomicsPreview.filteredRows().length,rows.filter(r=>["Up","Down"].includes(r.direction)).length);
  change("protein-filter","Not_tested");
  const notTested=w.proteomicsPreview.filteredRows();
  if(notTested.length){
    get("protein-rows").querySelector("button").click();await tick();
    assert.equal(get("protein-direction").textContent,"未检验");
  }
  change("protein-filter","all");
  const next=payload.samples.at(-1);
  click("nav-sample");await tick();assert.equal(get("sample-page").hidden,false);
  listeners.get("pca:plotly_click")({points:[{customdata:next.sample}]});await tick();
  assert.equal(current().sample,next.sample);assert.equal(get("sample-name").textContent,next.sample);
  assert.equal(get("sample-distribution")._traces[0].x.length,next.observed);
  change("sample-select",payload.samples[0].sample);await tick();assert.equal(current().sample,payload.samples[0].sample);
  get("sample-table").querySelectorAll("button")[1].click();await tick();assert.equal(current().sample,payload.samples[1].sample);
  click("nav-protein");await tick();assert.equal(get("protein-page").hidden,false);
  listeners.get("volcano:plotly_click")({points:[{customdata:target.protein_id}]});await tick();
  const abundanceSample=get("protein-abundance")._traces.find(t=>t.mode==="markers"&&t.customdata.length).customdata[0];
  listeners.get("protein-abundance:plotly_click")({points:[{customdata:abundanceSample}]});await tick();
  assert.equal(current().view,"sample");assert.equal(current().sample,abundanceSample);
  click("export-pca");click("export-volcano");
  assert.deepEqual(calls.map(c=>c.node),["pca","volcano"]);
  if(payload.contrasts.length>1){
    click("nav-protein");await tick();
    change("contrast-select",payload.contrasts[1].name);await tick();
    assert.equal(current().contrast,payload.contrasts[1].name);
    assert(get("comparison-note").textContent.includes("独立组"));
    assert(!get("protein-abundance")._traces.some(t=>t.mode==="lines"));
    change("contrast-select",co.name);await tick();
    const pairedTarget=rows.find(r=>r.n_complete_pairs===3);
    listeners.get("volcano:plotly_click")({points:[{customdata:pairedTarget.protein_id}]});await tick();
    assert.equal(get("protein-abundance")._traces.filter(t=>t.mode==="lines").length,3);
  }
  assert.equal(d.querySelectorAll("script[src],link[href]").length,0,"No network assets");
  await w.proteomicsPreview.setContrast(co.name);
  await w.proteomicsPreview.showView("qc");
  assert.equal(get("qc-bars")._traces[0].x.length,payload.samples.length);
  listeners.get("qc-bars:plotly_click")({points:[{customdata:next.sample}]});await tick();
  assert.equal(current().view,"sample");assert.equal(current().sample,next.sample);
  await w.proteomicsPreview.showView("heatmap");
  if(name==="DIA_cat")assert(get("protein-heatmap")._layout.annotations[0].text.includes("没有蛋白"));
  change("hm-set","variance");await tick();
  assert.equal(get("protein-heatmap")._traces[0].z.length,40);
  const heatProtein=get("protein-heatmap")._traces[0].customdata[0][0][0];
  listeners.get("protein-heatmap:plotly_click")({points:[{customdata:[heatProtein,1]}]});await tick();
  assert.equal(current().protein,heatProtein);assert.equal(current().view,"protein");
  await w.proteomicsPreview.showView("ora");change("ora-significance","all");await tick();
  const orows=w.proteomicsExtensions.oraRows();
  if(name==="DIA_cat")assert.equal(orows.length,0);
  if(orows.length){
    listeners.get("ora-plot:plotly_click")({points:[{customdata:orows[0].term}]});await tick();
    assert.equal(get("ora-title").textContent,orows[0].term);
    assert(get("ora-proteins").querySelector("button[data-open-protein]"));
  }
  await w.proteomicsPreview.showView("gsea");
  assert(get("gsea-curve")._traces.length===3,"Standard enrichment curve has 3 tracks");
  const sel=w.proteomicsExtensions.getSelection().gsea;
  const gsea=payload.extra.gsea.KEGG[co.name],grow=gsea.rows.find(r=>r.pathway===sel);
  const yy=get("gsea-curve")._traces[0].y;
  const es=Math.max(...yy)>-Math.min(...yy)?Math.max(...yy):Math.min(...yy);
  assert(Math.abs(es-grow.ES)<1e-7,"JavaScript GSEA curve matches fgsea ES");
  const hitgene=get("gsea-curve")._traces[1].customdata[0];
  listeners.get("gsea-curve:plotly_click")({points:[{customdata:hitgene}]});await tick();
  assert(get("gsea-gene-focus").textContent.includes(hitgene));
  const gp=get("gsea-leading").querySelector("[data-open-protein]");
  if(gp){const expected=gp.dataset.openProtein;gp.click();await tick();assert.equal(current().protein,expected);assert.equal(current().view,"protein");}
  await w.proteomicsPreview.showView("gsva");
  assert(get("gsva-heatmap")._traces[0].z.length>0);
  const pathway=get("gsva-heatmap")._traces[0].y[0];
  listeners.get("gsva-heatmap:plotly_click")({points:[{y:pathway}]});await tick();
  assert.equal(get("gsva-title").textContent,pathway);
  const scores=payload.extra.gsva.KEGG.scores[pathway];
  for(const tr of get("gsva-abundance")._traces.filter(t=>t.mode==="markers"))tr.customdata.forEach((s,i)=>assert.equal(tr.y[i],scores[s]));
  await w.proteomicsPreview.showView("ppi");
  if(name==="DIA_cat")assert(get("ppi-plot")._layout.annotations[0].text.includes("没有"));
  else if(payload.extra.ppi[co.name].nodes.length){
    const gene=payload.extra.ppi[co.name].nodes[0].gene;
    listeners.get("ppi-plot:plotly_click")({points:[{customdata:gene}]});await tick();
    assert.equal(get("ppi-gene").textContent,gene);
  }
  await w.proteomicsPreview.showView("files");
  assert(get("file-list").querySelectorAll("a").length>50);
  change("file-search","GSEA","input");assert([...get("file-list").querySelectorAll("a")].every(a=>a.textContent.includes("GSEA")));
  if(payload.contrasts.length>1){
    await w.proteomicsPreview.showView("gsea");
    const select=d.querySelector(".shared-contrast");select.value=payload.contrasts[1].name;
    select.dispatchEvent(new w.Event("change",{bubbles:true}));await tick();
    assert.equal(current().contrast,payload.contrasts[1].name);
    assert([...d.querySelectorAll(".shared-contrast")].every(e=>e.value===current().contrast));
  }
  const result={project:name,proteins:payload.proteins.length,samples:payload.samples.length,
    volcano_trace_proteins:tested.length,status:"PASS",
    validation:"9 modules, cross-view links and GSEA ES checked with mocked Plotly renderer; partial Plotly trace/layout schema validation; no real-browser rendering assertion"};
  dom.window.close();return result;
}
(async()=>{
  const results=[];
  for(const n of ["TMT_HepG2","DIA_cat"])results.push(await check(n));
  fs.writeFileSync(root+"/interaction_validation.json",JSON.stringify(results,null,2));
  console.log(JSON.stringify(results,null,2));
  schemaWindow.close();
})().catch(e=>{console.error(e);process.exitCode=1;});
