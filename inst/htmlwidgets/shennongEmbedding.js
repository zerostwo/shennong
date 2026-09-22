/* Shennong orthographic embedding viewer. No network, GPU or CDN dependency.
 * Projection and 96-bin transparent painter match R/plot_embedding.R.
 */
HTMLWidgets.widget({
  name: 'shennongEmbedding', type: 'output',
  factory: function(el, width, height) {
    let scene, camera, initial, canvas, ctx, stage, fields, output, frame = null;
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
      Object.keys(fields || {}).forEach(k => { fields[k].value = Number(camera[k].toFixed(5)); });
      if (output) output.value = cameraText();
    }
    function draw() {
      if (!scene || !canvas) return;
      const size = Math.max(1, Math.min(stage.clientWidth, stage.clientHeight));
      const dpr = window.devicePixelRatio || 1;
      if (canvas.width !== Math.round(size * dpr)) {
        canvas.width = canvas.height = Math.round(size * dpr);
        canvas.style.width = canvas.style.height = size + 'px';
      }
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.fillStyle = scene.control.background;
      ctx.fillRect(0, 0, size, size);
      const xy = p => [(p[0] + 1) * size / 2, (1 - p[1]) * size / 2];
      const points = scene.xyz.map((p, i) => ({p: project(p), i: i}));
      const faces = [];
      Object.values(scene.surfaces).forEach(surface => {
        const v = surface.vertices;
        for (let i = 0; i < v.length; i += 3) {
          const a = project(v[i]), b = project(v[i + 1]), c = project(v[i + 2]);
          const u = b.map((x, j) => (x - a[j]) / (j < 2 ? camera.zoom : 1));
          const w = c.map((x, j) => (x - a[j]) / (j < 2 ? camera.zoom : 1));
          const norm = [u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]];
          const rim = Math.pow(1 - Math.abs(norm[2]) / Math.max(Math.hypot(...norm), 1e-15), 3);
          const alpha = Math.min(.85, scene.control.surface_alpha + scene.control.glow * rim * .5);
          const light = scene.control.glow * rim * .55;
          const col = rgb(surface.color).map(x => Math.round(x + (255 - x) * light));
          faces.push({v: [a, b, c], z: (a[2] + b[2] + c[2])/3,
            fill: 'rgba(' + col.join(',') + ',' + alpha + ')'});
        }
      });
      let lo = Infinity, hi = -Infinity;
      points.forEach(p => {lo = Math.min(lo, p.p[2]); hi = Math.max(hi, p.p[2]);});
      faces.forEach(f => {lo = Math.min(lo, f.z); hi = Math.max(hi, f.z);});
      const bin = z => hi === lo ? 0 : Math.min(95, Math.floor(95 * (z - lo) / (hi - lo)));
      const pb = Array.from({length: 96}, () => []), fb = Array.from({length: 96}, () => []);
      points.forEach(p => pb[bin(p.p[2])].push(p)); faces.forEach(f => fb[bin(f.z)].push(f));
      // pt_size is millimetres; preview assumes a five-inch square panel.
      const radius = scene.pt_size / 127 * size * .375;
      function dot(p, scale, alpha) {
        const q = xy(p.p); ctx.beginPath();
        ctx.arc(q[0], q[1], Math.max(.15, radius * scale), 0, Math.PI * 2);
        ctx.fillStyle = rgba(scene.point_colors[p.i], alpha); ctx.fill();
      }
      for (let i = 0; i < 96; i++) {
        fb[i].sort((a,b) => a.z-b.z).forEach(f => {
          const v = f.v.map(xy); ctx.beginPath(); ctx.moveTo(...v[0]);
          ctx.lineTo(...v[1]); ctx.lineTo(...v[2]); ctx.closePath(); ctx.fillStyle = f.fill; ctx.fill();
        });
        pb[i].sort((a,b) => a.p[2]-b.p[2]);
        if (scene.control.glow > 0) pb[i].forEach(p => dot(p, 3, .08 * scene.control.glow));
        pb[i].forEach(p => dot(p, 1, scene.control.point_alpha));
      }
      ctx.fillStyle = '#ffffff'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = Math.max(10, scene.label_size / 127 * size) + 'px sans-serif';
      scene.labels.forEach((label, i) => { const p = xy(project(scene.label_xyz[i])); ctx.fillText(label, ...p); });
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
      frame = null; disposed = false; scene = x.scene; camera = JSON.parse(JSON.stringify(x.camera));
      initial = JSON.parse(JSON.stringify(camera)); rotate = scene.control.auto_rotate;
      el.replaceChildren();
      Object.assign(el.style, {display: 'flex', flexDirection: 'column', background: scene.control.background,
        color: 'white', font: '12px sans-serif', minHeight: '420px'});
      element('div', scene.title || '3D embedding', el).style.cssText = 'padding:10px 12px;font-size:15px;';
      const main = element('div', null, el); main.style.cssText = 'display:flex;flex:1;min-height:220px;overflow:hidden;';
      stage = element('div', null, main); stage.style.cssText = 'flex:1;min-width:0;display:flex;align-items:center;justify-content:center;overflow:hidden;';
      canvas = element('canvas', null, stage); canvas.setAttribute('aria-label', 'Rotate 3D embedding: drag; shift-drag to pan; scroll to zoom');
      canvas.style.touchAction = 'none'; ctx = canvas.getContext('2d');
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
      el.shennongCamera = () => JSON.parse(JSON.stringify(camera));
      el.shennongProject = p => project(p);
      sync(); schedule();
    }
    return {renderValue: render, resize: function() {schedule();}};
  }
});
