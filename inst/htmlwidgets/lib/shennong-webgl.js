/* Local WebGL renderer shared by interactive views and 600 dpi PDF captures.
 * One vertex per measured cell. No particle duplication or synthetic cells.
 * Copyright Shennong authors. MIT license (see package LICENSE).
 */
(function(global) {
  'use strict';
  function Renderer(canvas, scene) {
    const gl = canvas.getContext('webgl', {alpha: false, antialias: true,
      preserveDrawingBuffer: true, premultipliedAlpha: false});
    if (!gl) throw new Error('WebGL is unavailable. Enable WebGL or use a Chrome/Chromium installation with software rendering.');
    const resources = [], programs = [], surfaces = [], targets = [];
    let width = 0, height = 0;
    const color = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i+2), 16)/255);
    function program(vertex, fragment) {
      function compile(type, text) {
        const shader = gl.createShader(type); gl.shaderSource(shader, text); gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader));
        resources.push(['shader', shader]); return shader;
      }
      const p = gl.createProgram();
      gl.attachShader(p, compile(gl.VERTEX_SHADER, vertex)); gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fragment)); gl.linkProgram(p);
      if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
      programs.push(p); return p;
    }
    function buffer(data) {
      const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b);
      gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(data), gl.STATIC_DRAW); resources.push(['buffer', b]); return b;
    }
    function attr(p, name, b, size) {
      const at = gl.getAttribLocation(p, name); if (at < 0) return;
      gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.enableVertexAttribArray(at); gl.vertexAttribPointer(at, size, gl.FLOAT, false, 0, 0);
    }
    const uniform = (p, name) => gl.getUniformLocation(p, name);
    const normalized = rows => rows.flatMap(p => p.map((v,i) => (v-scene.center[i])/scene.radius));
    const common = 'attribute vec3 position; uniform mat3 rotation; uniform float zoom; uniform vec2 pan;';
    const projection = 'vec3 p=rotation*position; gl_Position=vec4(p.xy*zoom+pan,-p.z*.1,1.);';
    const surfaceProgram = program(common + `
      attribute vec3 normal; varying vec3 n;
      void main(){n=rotation*normal; ${projection}}`, `
      precision highp float; varying vec3 n; uniform vec3 tint;
      uniform float opacity; uniform float glow; uniform float glass;
      void main(){
        vec3 norm=normalize(n); float facing=abs(norm.z);
        float rim=pow(1.-facing,mix(2.8,1.5,glass));
        float diffuse=.35+.65*abs(dot(norm,normalize(vec3(-.4,.6,1.))));
        float alpha=opacity*(.10+.25*diffuse)+rim*(opacity+glow*mix(1.,.65,glass));
        vec3 lit=mix(tint*.65,tint, diffuse);
        lit=mix(lit,vec3(1.),rim*.10);
        gl_FragColor=vec4(lit,alpha*(gl_FrontFacing?1.:.45));
      }`);
    const pointProgram = program(common + `
      attribute vec3 tint; varying vec3 c; uniform float pointSize;
      void main(){c=tint; ${projection} gl_PointSize=pointSize;}`, `
      precision highp float; varying vec3 c; uniform float opacity; uniform float halo; uniform float whiten;
      void main(){
        float r=length(gl_PointCoord-vec2(.5))*2.; if(r>1.) discard;
        float a=mix(1.-smoothstep(.25,1.,r),exp(-r*r*4.),halo)*opacity;
        gl_FragColor=vec4(mix(c,vec3(1.),whiten*(1.-halo)),a);
      }`);
    const screenVertex = 'attribute vec2 position; varying vec2 uv; void main(){uv=position*.5+.5;gl_Position=vec4(position,0.,1.);}';
    const blurProgram = program(screenVertex, `
      precision highp float; varying vec2 uv; uniform sampler2D image;
      uniform vec2 stepSize; uniform float threshold;
      vec3 sampleAt(vec2 p){vec3 c=texture2D(image,p).rgb;return c*smoothstep(threshold,threshold+.3,max(c.r,max(c.g,c.b)));}
      void main(){vec3 c=sampleAt(uv)*.227027;
        c+=(sampleAt(uv+stepSize*1.384615)+sampleAt(uv-stepSize*1.384615))*.316216;
        c+=(sampleAt(uv+stepSize*3.230769)+sampleAt(uv-stepSize*3.230769))*.070270;
        gl_FragColor=vec4(c,1.);}`);
    const compositeProgram = program(screenVertex, `
      precision highp float; varying vec2 uv; uniform sampler2D base; uniform sampler2D bloom;
      uniform vec3 background; uniform float strength;
      void main(){vec3 light=texture2D(base,uv).rgb+texture2D(bloom,uv).rgb*strength;
        gl_FragColor=vec4(background+light*(1.-background),1.);}`);
    const quad = buffer([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]);
    const pointPositions = buffer(normalized(scene.xyz));
    const pointColors = buffer(scene.point_colors.flatMap(color));
    // Accumulate shared-vertex normals, then interpolate them per fragment.
    // The old flat face color caused polygon seams and an opaque sticker look.
    Object.values(scene.surfaces).forEach(s => {
      if (!s.vertices.length) return;
      const p = normalized(s.vertices), sums = new Map(), keys = [];
      for (let i=0;i<p.length;i+=3) keys.push(p.slice(i,i+3).map(v=>v.toFixed(6)).join(','));
      for (let i=0;i<p.length;i+=9) {
        const a=p.slice(i,i+3),b=p.slice(i+3,i+6),c=p.slice(i+6,i+9);
        const u=b.map((x,j)=>x-a[j]),v=c.map((x,j)=>x-a[j]);
        const n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
        for(let k=0;k<3;k++){const key=keys[i/3+k], old=sums.get(key)||[0,0,0]; sums.set(key,old.map((x,j)=>x+n[j]));}
      }
      const normals=keys.flatMap(key=>{const n=sums.get(key),len=Math.hypot(...n)||1;return n.map(x=>x/len);});
      surfaces.push({p:buffer(p),n:buffer(normals),count:p.length/3,color:color(s.color)});
    });
    function resize(w,h) {
      if(w===width && h===height) return;
      const max=gl.getParameter(gl.MAX_TEXTURE_SIZE);
      if(w>max || h>max) throw new Error('Requested 600 dpi raster exceeds this browser texture limit ('+max+' px). Reduce physical figure size.');
      width=w;height=h;canvas.width=w;canvas.height=h;
      targets.forEach(t=>{gl.deleteFramebuffer(t.f);gl.deleteTexture(t.t);});targets.length=0;
      for(let i=0;i<3;i++){
        const t=gl.createTexture();gl.bindTexture(gl.TEXTURE_2D,t);
        gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);
        gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,w,h,0,gl.RGBA,gl.UNSIGNED_BYTE,null);
        const f=gl.createFramebuffer();gl.bindFramebuffer(gl.FRAMEBUFFER,f);gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.COLOR_ATTACHMENT0,gl.TEXTURE_2D,t,0);
        if(gl.checkFramebufferStatus(gl.FRAMEBUFFER)!==gl.FRAMEBUFFER_COMPLETE) throw new Error('Unable to allocate the requested WebGL export size.');
        targets.push({f:f,t:t});
      }
    }
    function rotation(camera) {
      const a=camera.azimuth*Math.PI/180,e=camera.elevation*Math.PI/180,r=camera.roll*Math.PI/180;
      const x=[Math.cos(a),-Math.sin(a),0],y=[-Math.sin(e)*Math.sin(a),-Math.sin(e)*Math.cos(a),Math.cos(e)];
      const z=[Math.cos(e)*Math.sin(a),Math.cos(e)*Math.cos(a),Math.sin(e)];
      const row1=x.map((v,i)=>Math.cos(r)*v-Math.sin(r)*y[i]),row2=x.map((v,i)=>Math.sin(r)*v+Math.cos(r)*y[i]);
      return [row1[0],row2[0],z[0],row1[1],row2[1],z[1],row1[2],row2[2],z[2]];
    }
    function view(p,camera){gl.uniformMatrix3fv(uniform(p,'rotation'),false,rotation(camera));gl.uniform1f(uniform(p,'zoom'),camera.zoom);gl.uniform2fv(uniform(p,'pan'),camera.pan);}
    function screen(p,target){gl.bindFramebuffer(gl.FRAMEBUFFER,target);gl.useProgram(p);attr(p,'position',quad,2);}
    function texture(p,name,t,unit){gl.activeTexture(gl.TEXTURE0+unit);gl.bindTexture(gl.TEXTURE_2D,t);gl.uniform1i(uniform(p,name),unit);}
    this.render=function(camera,w,h,dpi){
      resize(Math.round(w),Math.round(h));gl.viewport(0,0,width,height);
      gl.disable(gl.DEPTH_TEST);gl.disable(gl.CULL_FACE);
      gl.bindFramebuffer(gl.FRAMEBUFFER,targets[0].f);gl.clearColor(0,0,0,1);gl.clear(gl.COLOR_BUFFER_BIT);
      gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE);
      gl.useProgram(surfaceProgram);view(surfaceProgram,camera);
      gl.uniform1f(uniform(surfaceProgram,'opacity'),scene.control.surface_alpha);
      gl.uniform1f(uniform(surfaceProgram,'glow'),scene.control.glow);
      gl.uniform1f(uniform(surfaceProgram,'glass'),scene.style==='glass'?1:0);
      surfaces.forEach(s=>{attr(surfaceProgram,'position',s.p,3);attr(surfaceProgram,'normal',s.n,3);gl.uniform3fv(uniform(surfaceProgram,'tint'),s.color);gl.drawArrays(gl.TRIANGLES,0,s.count);});
      gl.useProgram(pointProgram);view(pointProgram,camera);attr(pointProgram,'position',pointPositions,3);attr(pointProgram,'tint',pointColors,3);
      gl.uniform1f(uniform(pointProgram,'whiten'),scene.values===null?.45:0.);
      const size=scene.pt_size*dpi/25.4;
      gl.uniform1f(uniform(pointProgram,'halo'),1);gl.uniform1f(uniform(pointProgram,'opacity'),scene.control.point_alpha*scene.control.glow*.08);
      gl.uniform1f(uniform(pointProgram,'pointSize'),Math.max(1,size*4));gl.drawArrays(gl.POINTS,0,scene.xyz.length);
      gl.uniform1f(uniform(pointProgram,'halo'),0);gl.uniform1f(uniform(pointProgram,'opacity'),scene.control.point_alpha);
      gl.uniform1f(uniform(pointProgram,'pointSize'),Math.max(.6,size));gl.drawArrays(gl.POINTS,0,scene.xyz.length);
      gl.disable(gl.BLEND);
      screen(blurProgram,targets[1].f);texture(blurProgram,'image',targets[0].t,0);
      gl.uniform2f(uniform(blurProgram,'stepSize'),dpi/100/width,0);gl.uniform1f(uniform(blurProgram,'threshold'),.16);gl.drawArrays(gl.TRIANGLES,0,6);
      screen(blurProgram,targets[2].f);texture(blurProgram,'image',targets[1].t,0);
      gl.uniform2f(uniform(blurProgram,'stepSize'),0,dpi/100/height);gl.uniform1f(uniform(blurProgram,'threshold'),0);gl.drawArrays(gl.TRIANGLES,0,6);
      screen(compositeProgram,null);texture(compositeProgram,'base',targets[0].t,0);texture(compositeProgram,'bloom',targets[2].t,1);
      gl.uniform3fv(uniform(compositeProgram,'background'),color(scene.control.background));gl.uniform1f(uniform(compositeProgram,'strength'),scene.control.glow*.65);
      gl.drawArrays(gl.TRIANGLES,0,6);
      const error=gl.getError();if(error!==gl.NO_ERROR) throw new Error('WebGL rendering failed: '+error);
    };
    this.capture=function(camera,w,h,dpi){this.render(camera,w,h,dpi);return canvas.toDataURL('image/png');};
    this.dispose=function(){targets.forEach(t=>{gl.deleteFramebuffer(t.f);gl.deleteTexture(t.t);});resources.forEach(r=>{if(r[0]==='shader')gl.deleteShader(r[1]);else gl.deleteBuffer(r[1]);});programs.forEach(p=>gl.deleteProgram(p));const ext=gl.getExtension('WEBGL_lose_context');if(ext)ext.loseContext();};
  }
  global.ShennongEmbeddingRenderer=Renderer;
})(window);
