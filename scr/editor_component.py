import streamlit as st

# Bidirectional component: 3Dmol.js handles atom picking, while Streamlit receives
# selected atom indices. No Node/npm build is needed with Streamlit components v2.
EDITOR_JS = r"""
export default function({data, parentElement, setStateValue, state}) {
  const root = parentElement.querySelector('.molecule-editor');
  const status = parentElement.querySelector('.selection-status');
  if (!root) return;
  // Remove the old viewer when Streamlit changes the molecule.
  // Save the current camera before rebuilding the WebGL canvas.
  if (root._gaussianViewer) {
    try { sessionStorage.setItem(data.cameraKey, JSON.stringify(root._gaussianViewer.getView())); }
    catch(e) {}
  }
  root.replaceChildren();
  const picked = Array.isArray(state?.selected) ? [...state.selected] : [];
  const maxPick = data.maxPick || 4;
  function display() {
    status.textContent = picked.length
      ? 'Selected atoms: ' + picked.map(i => (i + 1)).join(' → ')
      : 'Click atoms in the order requested by the editing operation.';
  }
  display();
  function build() {
    if (!window.$3Dmol) {
      status.textContent = 'Unable to load 3Dmol.js. Check internet access.';
      return;
    }
    const viewer = $3Dmol.createViewer(root, {backgroundColor: data.bg, antialias:true});
    const atoms = data.atoms.map((a,i) => ({
      elem:a.elem, x:a.x, y:a.y, z:a.z,
      serial:i+1, bonds:[], bondOrder:[]
    }));
    const hbonds=[];
    for(const [i,j] of data.bonds) {
      const kind=(data.bondTypes || {})[[Math.min(i,j),Math.max(i,j)].join('-')] || 'Single';
      if(kind==='H-bond'){hbonds.push([i,j]);continue;}
      const order=kind==='Double'?2:kind==='Triple'?3:kind==='Aromatic'?1.5:1;
      atoms[i].bonds.push(j); atoms[i].bondOrder.push(order);
      atoms[j].bonds.push(i); atoms[j].bondOrder.push(order);
    }
    root._gaussianViewer = viewer;
    const model = viewer.addModel();
    model.addAtoms(atoms);
    hbonds.forEach(([i,j])=>viewer.addLine({start:{x:atoms[i].x,y:atoms[i].y,z:atoms[i].z},
      end:{x:atoms[j].x,y:atoms[j].y,z:atoms[j].z},dashed:true,color:'#22d3ee',linewidth:2}));
    viewer.setStyle({}, {stick:{radius:0.15}, sphere:{scale:0.27}});
    // Preview internal-coordinate edits against the original committed geometry.
    // The Python edit_geometry function remains authoritative when Apply is clicked.
    function preview() {
      // Rebuild from committed coordinates to avoid cumulative slider drift.
      const base=model.selectedAtoms({});
      base.forEach((a,i)=>{a.x=data.atoms[i].x;a.y=data.atoms[i].y;a.z=data.atoms[i].z;});
      atoms.forEach((a,i)=>{a.x=data.atoms[i].x;a.y=data.atoms[i].y;a.z=data.atoms[i].z;});
      const ids=picked;
      const op=data.operation;
      const val=Number(data.value);
      if (!Number.isFinite(val) || ids.length !== data.maxPick ||
          !['Bond length','Bond angle','Dihedral angle'].includes(op)) return;
      const pts=data.atoms.map(a=>[a.x,a.y,a.z]);
      const sub=(a,b)=>a.map((x,i)=>x-b[i]);
      const add=(a,b)=>a.map((x,i)=>x+b[i]);
      const mul=(a,t)=>a.map(x=>x*t);
      const dot=(a,b)=>a.reduce((t,x,i)=>t+x*b[i],0);
      const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
      const norm=a=>Math.sqrt(dot(a,a));
      const unit=a=>{const n=norm(a);if(n<1e-9)throw Error('Degenerate axis');return mul(a,1/n)};
      const edges=data.bonds;
      const side=(start,cutA,cutB)=>{
        const seen=new Set([start]), queue=[start];
        while(queue.length){const i=queue.pop();for(const [a,b] of edges){
          if((a===cutA&&b===cutB)||(a===cutB&&b===cutA))continue;
          const j=a===i?b:b===i?a:-1;
          if(j>=0&&!seen.has(j)){seen.add(j);queue.push(j)}
        }}return seen;
      };
      const rotate=(set,origin,axis,deg)=>{
        const k=unit(axis),t=deg*Math.PI/180;
        for(const i of set){const v=sub(pts[i],origin);
          pts[i]=add(origin,add(add(mul(v,Math.cos(t)),mul(cross(k,v),Math.sin(t))),mul(k,dot(k,v)*(1-Math.cos(t)))));}
      };
      try {
        if(op==='Bond length'){
          const [a,b]=ids, group=side(b,a,b);
          if(group.has(a))return;
          const vec=sub(pts[b],pts[a]);
          const shift=mul(unit(vec),val-norm(vec));
          for(const i of group)pts[i]=add(pts[i],shift);
        } else if(op==='Bond angle'){
          const [a,b,c]=ids,group=side(c,b,c);
          if(group.has(a)||group.has(b))return;
          const u=unit(sub(pts[a],pts[b])),v=unit(sub(pts[c],pts[b]));
          const current=Math.acos(Math.max(-1,Math.min(1,dot(u,v))))*180/Math.PI;
          let axis=cross(v,u);
          if(norm(axis)<1e-9)axis=cross(v,[1,0,0]);
          if(norm(axis)<1e-9)axis=cross(v,[0,1,0]);
          rotate(group,pts[b],axis,current-val);
        } else {
          const [a,b,c,d]=ids,group=side(d,b,c);
          if(group.has(a)||group.has(b))return;
          const axis=unit(sub(pts[c],pts[b]));
          const proj=v=>sub(v,mul(axis,dot(v,axis)));
          const v=unit(proj(sub(pts[a],pts[b]))),w=unit(proj(sub(pts[d],pts[c])));
          const current=Math.atan2(dot(cross(v,w),axis),dot(v,w))*180/Math.PI;
          rotate(group,pts[c],axis,val-current);
        }
        const modelAtoms=model.selectedAtoms({});
        modelAtoms.forEach((a,i)=>{a.x=pts[i][0];a.y=pts[i][1];a.z=pts[i][2]});
        model.setCoordinates(modelAtoms.map(a=>[a.x,a.y,a.z]),'array');
        atoms.forEach((a,i)=>{a.x=pts[i][0];a.y=pts[i][1];a.z=pts[i][2]});
      } catch(e) {console.warn('Preview unavailable:',e)}
    }
    // Geometry changes are committed by Python and displayed on the next rerun.
    function highlight() {
      viewer.removeAllLabels();
      viewer.setStyle({}, {stick:{radius:0.15}, sphere:{scale:0.27}});
      picked.forEach((idx, position) => {
        viewer.addStyle({serial:idx+1}, {sphere:{radius:0.43, color:'#facc15'}});
        viewer.addLabel(String(position+1), {
          position:{x:atoms[idx].x, y:atoms[idx].y+0.55, z:atoms[idx].z},
          fontColor:'#111827', backgroundColor:'#fde68a', backgroundOpacity:0.95,
          fontSize:15, borderThickness:0, inFront:true
        });
      });
      viewer.render();
    }
    viewer.setClickable({}, true, atom => {
      const idx = atom.serial - 1;
      const old = picked.indexOf(idx);
      if(old >= 0) picked.splice(old,1);
      else {
        if(picked.length >= maxPick) picked.shift();
        picked.push(idx);
      }
      display(); highlight();
      setStateValue('selected', picked.slice());
    });
    viewer.zoomTo();
    try {
      const saved = JSON.parse(sessionStorage.getItem(data.cameraKey) || 'null');
      if (saved) viewer.setView(saved);
    } catch(e) {}
    const saveCamera = () => {
      try { sessionStorage.setItem(data.cameraKey, JSON.stringify(viewer.getView())); }
      catch(e) {}
    };
    root.addEventListener('pointerup', saveCamera);
    root.addEventListener('wheel', saveCamera);
    root.addEventListener('touchend', saveCamera);
    highlight(); display();
    const ro = new ResizeObserver(() => viewer.resize());
    ro.observe(root);
    return () => { ro.disconnect(); viewer.clear(); };
  }
  if(window.$3Dmol) return build();
  const script = document.createElement('script');
  script.src = 'https://3Dmol.org/build/3Dmol-min.js';
  script.onload = build;
  script.onerror = () => {status.textContent='3Dmol.js could not load.'};
  document.head.appendChild(script);
}
"""

PICKER = st.components.v2.component(
    'gaussian_explorer_atom_picker',
    html='<div class="molecule-editor"></div><div class="selection-status"></div>',
    css=""".molecule-editor {width:100%;height:530px;border-radius:14px;overflow:hidden;}
.selection-status {padding:9px 12px;font-size:14px;color:var(--st-text-color);}
""",
    js=EDITOR_JS,
)


