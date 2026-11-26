import streamlit as st
import streamlit.components.v1 as components

def viewer(xyz=None, frames=None, cube=None, style='Ball and stick', bg='#AABDE1',
           height=560, spin=False, animate=False, amplitude=1.0, iso=0.035):
    """Embed 3Dmol.js using a controlled HTML document. No user text is injected as JS."""
    import json
    models = frames if frames is not None else [xyz]
    models = [m for m in models if m]
    if not models:
        st.info('No molecular coordinates available.')
        return
    # Escape JSON in HTML script context.
    payload = json.dumps({'models': models, 'cube': cube, 'style': style,
                          'bg': bg, 'spin': spin, 'animate': animate, 'iso': iso}).replace('</', '<\\/')
    html = '''<!doctype html><html><head><script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>html,body,#viewer{height:100%;width:100%;margin:0;background:#AABDE1;overflow:hidden}</style></head>
<body><div id="viewer"></div><script>
const p=PAYLOAD;
function start(){
 const v=$3Dmol.createViewer('viewer',{backgroundColor:p.bg,antialias:true});
 if(p.models.length>1){v.addModelsAsFrames(p.models.join(''),'xyz');}
 else {v.addModel(p.models[0],'xyz');}
 const sty=p.style==='Stick'?{stick:{radius:0.16}}:p.style==='Sphere'?{sphere:{scale:0.32}}:
   p.style==='Line'?{line:{linewidth:2}}:{stick:{radius:0.14},sphere:{scale:0.27}};
 v.setStyle({},sty);
 if(p.cube){
   v.addVolumetricData(p.cube,'cube',{isoval:p.iso,color:'#38bdf8',opacity:0.75});
   v.addVolumetricData(p.cube,'cube',{isoval:-p.iso,color:'#fb7185',opacity:0.75});
 }
 v.zoomTo(); if(p.spin)v.spin('y',0.4);
 if(p.animate && p.models.length>1)v.animate({loop:'backAndForth',reps:0,interval:110});
 v.render();
}
if(window.$3Dmol){start()}else{document.body.innerHTML='3Dmol.js could not load. Check internet access.'}
</script></body></html>'''.replace('PAYLOAD',payload)
    components.html(html, height=height, scrolling=False)
