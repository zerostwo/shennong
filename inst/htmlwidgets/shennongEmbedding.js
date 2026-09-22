/* Local WebGL embedding viewer. The same renderer supplies static PDF rasters. */
HTMLWidgets.widget({
  name: 'shennongEmbedding', type: 'output',
  factory: function(el, width, height) {
    let scene, camera, initial, canvas, labelCanvas, ctx, renderer, stage, fields, output, frame = null;
    let rotate = false, lastTime = 0, disposed = false;
    const rad = Math.PI / 180;
    function project(p) {
      const q = p.map((v, i) => (v - scene.center[i]) / scene.radius);
      const a = camera.azimuth * rad, e = camera.elevation * rad, r = camera.roll * rad;
      const x = Math.cos(a) * q[0] - Math.sin(a) * q[1];
      const along = Math.sin(a) * q[0] + Math.cos(a) * q[1];
      const y = -Math.sin(e) * along + Math.cos(e) * q[2];
      return [(Math.cos(r) * x - Math.sin(r) * y) * camera.zoom + camera.pan[0],
        (Math.sin(r) * x + Math.cos(r) * y) * camera.zoom + camera.pan[1],
        Math.cos(e) * along + Math.sin(e) * q[2]];
    }
    function rgb(hex) {
      const m = /^#([0-9a-f]{6})/i.exec(hex);
      if (!m) return [180, 180, 180];
      return [0, 2, 4].map(i => parseInt(m[1].slice(i, i + 2), 16));
    }
    function rgba(color, alpha) { return 'rgba(' + rgb(color).join(',') + ',' + alpha + ')'; }
    function cameraText() {
      const n = x => Number(x.toFixed(8));
      return 'camera <- list(azimuth = ' + n(camera.azimuth) + ', elevation = ' + n(camera.elevation) +
        ', roll = ' + n(camera.roll) + ', zoom = ' + n(camera.zoom) +
        ', pan = c(' + camera.pan.map(n).join(', ') + '))';
    }
    function sync() {
      Object.keys(fields || {}).forEach(k => {
        // A pending GPU frame must not overwrite an angle the user is typing.
        if (document.activeElement !== fields[k]) fields[k].value = Number(camera[k].toFixed(5));
      });
      if (output) output.value = cameraText();
    }
    function draw() {
      if (!scene || !renderer) return;
      const size = Math.max(1, Math.min(stage.clientWidth, stage.clientHeight));
      const dpr = window.devicePixelRatio || 1;
      const pixels = Math.round(size * dpr);
      canvas.style.width = canvas.style.height = size + 'px';
      labelCanvas.width = labelCanvas.height = pixels;
      labelCanvas.style.width = labelCanvas.style.height = size + 'px';
      renderer.render(camera, pixels, pixels, pixels / 5);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, size, size);
      const textSize = scene.label_size * size / 127;
      ctx.font = Math.max(10, textSize) + 'px sans-serif';
      ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
      scene.labels.forEach((label, i) => {
        const q = project(scene.label_xyz[i]), x = (q[0]+1)*size/2, y = (1-q[1])*size/2;
        const tint = scene.group_color_values[scene.group_labels.indexOf(label)];
        const w = ctx.measureText(label).width;
        ctx.fillStyle = 'rgba(6,9,20,.72)';ctx.fillRect(x+5,y-9,w+10,18);
        ctx.fillStyle = tint;
        if(scene.style==='nebula') ctx.fillRect(x,y-9,2,18);
        else {ctx.shadowColor=tint;ctx.shadowBlur=6;ctx.beginPath();ctx.arc(x,y,3,0,Math.PI*2);ctx.fill();ctx.shadowBlur=0;}
        ctx.fillStyle = '#f1f4fa';ctx.fillText(label,x+10,y);
      });
      sync();
    }
    function schedule() {
      if (frame !== null) return;
      frame = requestAnimationFrame(function tick(time) {
        frame = null;
        if (disposed || !el.isConnected) return;
        if (rotate) {camera.azimuth += lastTime ? Math.min(time-lastTime, 100) * .008 : 0; lastTime = time;}
        else lastTime = 0;
        draw(); if (rotate) schedule();
      });
    }
    function element(tag, text, parent) {
      const node = document.createElement(tag); if (text !== null) node.textContent = text;
      if (parent) parent.appendChild(node); return node;
    }
    function pause() { rotate = false; if (el.querySelector('[data-rotate]')) el.querySelector('[data-rotate]').checked = false; }
    function render(x) {
      if (frame !== null) cancelAnimationFrame(frame);
      if (renderer) renderer.dispose();
      frame = null; disposed = false; scene = x.scene; camera = JSON.parse(JSON.stringify(x.camera));
      initial = JSON.parse(JSON.stringify(camera)); rotate = scene.control.auto_rotate;
      el.replaceChildren();
      Object.assign(el.style, {display: 'flex', flexDirection: 'column', background: scene.control.background,
        color: 'white', font: '12px sans-serif', minHeight: '420px'});
      element('div', scene.title || '3D embedding', el).style.cssText = 'padding:10px 12px;font-size:15px;';
      const main = element('div', null, el); main.style.cssText = 'display:flex;flex:1;min-height:220px;overflow:hidden;';
      stage = element('div', null, main); stage.style.cssText = 'flex:1;min-width:0;display:flex;align-items:center;justify-content:center;overflow:hidden;';
      const canvasBox = element('div', null, stage); canvasBox.style.cssText = 'position:relative;display:flex;';
      canvas = element('canvas', null, canvasBox); canvas.setAttribute('aria-label', 'Rotate 3D embedding: drag; shift-drag to pan; scroll to zoom');
      canvas.style.touchAction = 'none';
      labelCanvas = element('canvas', null, canvasBox); labelCanvas.style.cssText = 'position:absolute;left:0;top:0;pointer-events:none;';
      ctx = labelCanvas.getContext('2d');
      try {
        renderer = new ShennongEmbeddingRenderer(canvas, scene);
      } catch (error) {
        el.shennongError = error.message;
        const alert = element('div', error.message, el); alert.setAttribute('role', 'alert');
        throw error;
      }
      if (scene.show_legend) {
        const legend = element('div', null, main); legend.style.cssText = 'width:145px;flex-shrink:0;align-self:center;padding:8px;max-height:100%;overflow:auto;';
        element('strong', scene.legend_title, legend);
        if (scene.values === null) {
          scene.group_labels.forEach((name, i) => {
            const row = element('div', null, legend); row.style.paddingTop = '6px';
            const dot = element('span', '● ', row); dot.style.color = scene.group_color_values[i]; element('span', name, row);
          });
        } else {
          const gradient = element('div', null, legend);
          gradient.style.cssText = 'height:12px;margin:8px 0;background:linear-gradient(to right,' + scene.ramp.join(',') + ');';
          element('div', scene.limits.map(v => Number(v.toPrecision(4))).join(' — '), legend);
        }
      }
      const controls = element('div', null, el); controls.style.cssText = 'padding:8px 12px;border-top:1px solid #303748;';
      element('div', 'Drag to rotate · Shift-drag to pan · Scroll to zoom', controls).style.marginBottom = '6px';
      fields = {};
      ['azimuth', 'elevation', 'roll', 'zoom'].forEach(k => {
        const label = element('label', k + ' ', controls); label.style.marginRight = '8px';
        const input = element('input', null, label); input.type = 'number'; input.step = k === 'zoom' ? '.05' : '1';
        input.style.width = '70px'; input.setAttribute('aria-label', k); fields[k] = input;
        input.onfocus = pause;
        input.onchange = () => {const v = Number(input.value); if (Number.isFinite(v) && (k !== 'zoom' || (v >= .01 && v <= 100))) {pause(); camera[k] = v; schedule();}};
      });
      const rotation = element('label', ' Auto rotate ', controls), checkbox = element('input', null, rotation);
      checkbox.type = 'checkbox'; checkbox.dataset.rotate = 'true'; checkbox.checked = rotate;
      checkbox.onchange = () => {rotate = checkbox.checked; schedule();};
      const buttons = element('div', null, controls); buttons.style.marginTop = '8px';
      function button(name, fn) {const b = element('button', name, buttons); b.style.marginRight = '8px'; b.onclick = fn;}
      button('Reset view', () => {pause(); camera = JSON.parse(JSON.stringify(initial)); schedule();});
      button('Copy R camera', async () => {
        pause(); sync(); output.select();
        try {await navigator.clipboard.writeText(output.value);} catch (e) {document.execCommand('copy');}
      });
      button('Download camera JSON', () => {
        pause(); sync();
        const url = URL.createObjectURL(new Blob([JSON.stringify(camera, null, 2)], {type: 'application/json'}));
        const a = element('a', null, el); a.href = url; a.download = 'shennong-camera.json'; a.click(); a.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
      });
      output = element('textarea', null, controls); output.readOnly = true;
      output.setAttribute('aria-label', 'R camera parameters'); output.style.cssText = 'display:block;width:98%;height:36px;margin-top:8px;';
      let drag = null;
      canvas.onpointerdown = event => {pause(); drag = {x:event.clientX,y:event.clientY}; canvas.setPointerCapture(event.pointerId);};
      canvas.onpointermove = event => {
        if (!drag) return;
        const dx = event.clientX-drag.x, dy = event.clientY-drag.y;
        if (event.shiftKey) {camera.pan[0] += 2*dx/canvas.clientWidth; camera.pan[1] -= 2*dy/canvas.clientHeight;}
        else {camera.azimuth += dx*.4; camera.elevation -= dy*.4;}
        drag = {x:event.clientX,y:event.clientY}; schedule();
      };
      canvas.onpointerup = canvas.onpointercancel = () => {drag = null;};
      canvas.onwheel = event => {event.preventDefault(); pause(); camera.zoom = Math.max(.01, Math.min(100, camera.zoom * Math.exp(-event.deltaY*.001))); schedule();};
      el.shennongCapture = (w,h,dpi) => renderer.capture(camera,w,h,dpi);
      el.shennongCamera = () => JSON.parse(JSON.stringify(camera));
      el.shennongProject = p => project(p);
      sync(); schedule();
    }
    return {renderValue: render, resize: function() {schedule();}};
  }
});
