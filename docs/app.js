/* Norilsk Ground Truth: renders docs/data/site.json (built by 13_build_site.py). */
const ZONES = ["industrial", "accident", "urban", "control_NW", "control_SE", "reference"];
const ZLABEL = {
  industrial: "Industrial", accident: "TPP-3 site", urban: "City",
  control_NW: "Tundra NW", control_SE: "Tundra SE", reference: "Plateau bedrock"
};
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const ZCOL = Object.fromEntries(ZONES.map((z) => [z, css(`--z-${z}`)]));
const fmt = (v, d = 1) => (v == null || !isFinite(v) ? "–" : (+v).toFixed(d));
const tip = document.getElementById("tooltip");

function showTip(html, ev) {
  tip.innerHTML = html;
  tip.style.opacity = 1;
  const pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
  let x = ev.clientX + pad, y = ev.clientY + pad;
  if (x + w > innerWidth - 8) x = ev.clientX - w - pad;
  if (y + h > innerHeight - 8) y = ev.clientY - h - pad;
  tip.style.left = x + "px"; tip.style.top = y + "px";
}
const hideTip = () => (tip.style.opacity = 0);

fetch("data/site.json").then((r) => r.json()).then(init).catch((e) => {
  console.error(e);
  document.getElementById("hero-stats").innerHTML =
    `<div class="stat"><div class="v">–</div><div class="k">Data not built yet: run 13_build_site.py</div></div>`;
});

function init(S) {
  window.SITE = S;
  // each section is independent: one failure (e.g. no WebGL) must not blank the rest
  const safe = (fn, name) => { try { fn(S); } catch (e) { console.error(name, e); return e; } };
  safe(heroStats, "hero"); safe(timeline, "timeline");
  if (safe(mapInit, "map")) document.getElementById("map").innerHTML =
    `<p class="note" style="padding:24px">The interactive map needs WebGL, which this browser does not provide. The figures below still work.</p>`;
  safe(results, "results"); safe(local, "local"); safe(verdict, "verdict"); safe(nisar, "nisar");
  document.getElementById("generated").textContent = `Data build: ${S.generated}`;
  reveal();
}

/* ---------- hero ---------- */
function heroStats(S) {
  const n1 = Object.values(S.timeline).reduce((a, t) => a + t.s1_dates, 0);
  const proc = Object.values(S.s1_meta || {}).reduce((a, m) => a + m.pairs.length, 0);
  const tiles = [
    { v: n1, k: "Sentinel-1 summer passes over Norilsk, 2017–2026" },
    { v: proc || "30", k: "interferograms inverted for summers 2020 and 2025" },
    { v: Object.keys((S.nisar && S.nisar.products) || {}).length || 4, k: "NISAR L-band pairs, first science summer 2026" }
  ];
  const h = headline(S);
  if (h) tiles.push({ v: h.v, small: h.small, k: h.k, hl: true });
  document.getElementById("hero-stats").innerHTML = tiles.map((t) =>
    `<div class="stat ${t.hl ? "hl" : ""}"><div class="v">${t.v}${t.small ? `<small>${t.small}</small>` : ""}</div><div class="k">${t.k}</div></div>`).join("");
}

function fitOf(S, zone, yr) { return S.fits.find((f) => f.zone === zone && +f.summer === +yr); }

function headline(S) {
  const L = S.local;
  if (L) return { v: fmt(L.r_local, 2), small: "r", k: `correlation of local thaw subsidence between 2020 and 2025 (chance: ${fmt(L.r_null_shift_mean, 2)} ± ${fmt(L.r_null_shift_std, 2)})` };
  const V = S.verdict;
  if (!V || V.status !== "ready") return null;
  const n = V.n_summers, k = V.anthropogenic.yes;
  return { v: `${k}/${n}`, k: "summers in which industrial ground sank faster than natural tundra beyond 2σ" };
}

/* ---------- timeline ---------- */
function timeline(S) {
  const el = document.getElementById("timeline");
  const years = Object.keys(S.timeline).map(Number);
  const W = 1100, H = 170, m = { l: 10, r: 10, t: 24, b: 44 };
  const x = d3.scaleBand().domain(years).range([m.l, W - m.r]).padding(0.18);
  const maxN = d3.max(years, (y) => S.timeline[y].s1_dates) || 10;
  const y = d3.scaleLinear().domain([0, maxN]).range([H - m.b, m.t + 18]);
  const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
  const defs = svg.append("defs");
  defs.append("pattern").attr("id", "hatch").attr("width", 8).attr("height", 8)
    .attr("patternUnits", "userSpaceOnUse").attr("patternTransform", "rotate(45)")
    .append("line").attr("x1", 0).attr("y1", 0).attr("x2", 0).attr("y2", 8)
    .attr("stroke", "#3a3f48").attr("stroke-width", 3);
  svg.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", H - m.b).attr("y2", H - m.b).attr("stroke", "#383835");
  years.forEach((yr) => {
    const t = S.timeline[yr], g = svg.append("g");
    const bx = x(yr), bw = x.bandwidth();
    if (t.s1_dates === 0 && !t.nisar) {
      g.append("rect").attr("x", bx).attr("y", m.t + 18).attr("width", bw).attr("height", H - m.b - m.t - 18)
        .attr("fill", "url(#hatch)").attr("rx", 4);
    }
    if (t.s1_dates > 0) {
      g.append("rect").attr("x", bx + bw * 0.12).attr("width", bw * 0.76).attr("y", y(t.s1_dates))
        .attr("height", H - m.b - y(t.s1_dates)).attr("rx", 4)
        .attr("fill", t.processed ? css("--z-industrial") : "#3a3f48");
      g.append("text").attr("x", bx + bw / 2).attr("y", y(t.s1_dates) - 6).attr("text-anchor", "middle")
        .attr("fill", t.processed ? css("--ink") : css("--muted")).attr("font-size", 12).text(t.s1_dates);
    }
    if (t.nisar) {
      g.append("rect").attr("x", bx + bw * 0.12).attr("width", bw * 0.76).attr("y", m.t).attr("height", 12).attr("rx", 3)
        .attr("fill", css("--accent"));
    }
    g.append("text").attr("x", bx + bw / 2).attr("y", H - m.b + 20).attr("text-anchor", "middle")
      .attr("fill", css("--ink-2")).attr("font-size", 13).attr("font-weight", 600).text(yr);
    g.append("rect").attr("x", bx).attr("y", 0).attr("width", bw).attr("height", H).attr("fill", "transparent")
      .on("mousemove", (ev) => showTip(`<b>${yr}</b><div class="tt-row"><span class="tt-k">Sentinel-1 summer passes</span><span>${t.s1_dates}</span></div>` +
        (t.processed ? `<div class="tt-k">processed in this study</div>` : "") +
        (t.nisar ? `<div class="tt-k">NISAR L-band, both directions</div>` : "") +
        (t.s1_dates === 0 && !t.nisar ? `<div class="tt-k">no radar acquisitions</div>` : ""), ev))
      .on("mouseleave", hideTip);
  });
  const gap = years.filter((yy) => S.timeline[yy].s1_dates === 0 && !S.timeline[yy].nisar);
  if (gap.length) {
    const gx = (x(gap[0]) + x(gap[gap.length - 1]) + x.bandwidth()) / 2;
    svg.append("text").attr("x", gx).attr("y", (H - m.b + m.t + 18) / 2 + 4).attr("text-anchor", "middle")
      .attr("fill", css("--ink-2")).attr("font-size", 13).text("no Sentinel-1 data");
  }
  svg.append("text").attr("x", W - m.r).attr("y", H - 6).attr("text-anchor", "end").attr("fill", css("--muted")).attr("font-size", 12)
    .text("bars: Sentinel-1 summer passes (orange = processed) · gold strip: NISAR");
  document.getElementById("timeline-note").textContent =
    "Sentinel-1B failed in December 2021 and Sentinel-1A did not cover Norilsk afterwards; Sentinel-1C restored coverage in 2025. NISAR's provisional products start in June 2026.";
}

/* ---------- map ---------- */
let MAP, CURRENT;
function mapInit(S) {
  const style = {
    version: 8,
    glyphs: "https://demotiles.maplibre.org/font/{fontstack}/{range}.pbf",
    sources: {
      sat: { type: "raster", tileSize: 256, maxzoom: 17,
        tiles: ["https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"],
        attribution: "Imagery © Esri, Maxar, Earthstar Geographics" },
      dark: { type: "raster", tileSize: 256,
        tiles: ["https://a.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png", "https://b.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png"],
        attribution: "© CARTO © OpenStreetMap contributors" }
    },
    layers: [
      { id: "sat", type: "raster", source: "sat", paint: { "raster-saturation": -0.35, "raster-brightness-max": 0.7 } },
      { id: "dark", type: "raster", source: "dark", layout: { visibility: "none" } }
    ]
  };
  MAP = new maplibregl.Map({ container: "map", style, center: [88.2, 69.39], zoom: 9.3, attributionControl: { compact: true },
    cooperativeGestures: matchMedia("(max-width: 900px)").matches });
  MAP.addControl(new maplibregl.NavigationControl({ showCompass: false }), "top-left");
  MAP.addControl(new maplibregl.ScaleControl({ unit: "metric" }), "bottom-left");

  const keys = Object.keys(S.overlays);
  const order = ["local_persistent", "local_2020", "local_2025", "s1_2020", "s1_2025", "nisar_late_vertical", "nisar_late_east", "nisar_early_vertical", "nisar_early_east"];
  const ordered = order.filter((k) => keys.includes(k)).concat(keys.filter((k) => !order.includes(k)));
  const short = { local_persistent: "Repeats ★", local_2020: "Local 2020", local_2025: "Local 2025", s1_2020: "S1 2020", s1_2025: "S1 2025", nisar_early_vertical: "NISAR early ↕", nisar_early_east: "NISAR early ↔",
    nisar_late_vertical: "NISAR late ↕", nisar_late_east: "NISAR late ↔" };
  const seg = document.getElementById("layer-seg");
  seg.innerHTML = ordered.map((k, i) => `<button role="radio" aria-checked="${i === 0}" data-k="${k}">${short[k] || k}</button>`).join("");

  MAP.on("load", () => {
    ordered.forEach((k, i) => {
      const o = S.overlays[k];
      MAP.addSource(k, { type: "image", url: o.file, coordinates: o.coords });
      MAP.addLayer({ id: k, type: "raster", source: k, paint: { "raster-opacity": 0.92, "raster-fade-duration": 250, "raster-resampling": "linear" },
        layout: { visibility: i === 0 ? "visible" : "none" } });
    });
    fetch("data/zones.geojson").then((r) => r.json()).then((gj) => {
      MAP.addSource("zones", { type: "geojson", data: gj });
      const colorExpr = ["match", ["get", "zone"], ...ZONES.flatMap((z) => [z, ZCOL[z]]), "#888"];
      MAP.addLayer({ id: "zones-fill", type: "fill", source: "zones", paint: { "fill-color": colorExpr, "fill-opacity": 0.04 } });
      MAP.addLayer({ id: "zones-line", type: "line", source: "zones", paint: { "line-color": colorExpr, "line-width": 1.6, "line-opacity": 0.95 } });
      MAP.addLayer({ id: "zones-hl", type: "line", source: "zones", paint: { "line-color": "#ffffff", "line-width": 2.6 }, filter: ["==", "zone", ""] });
      MAP.on("click", "zones-fill", (e) => selectZone(e.features[0].properties.zone));
      MAP.on("mouseenter", "zones-fill", () => (MAP.getCanvas().style.cursor = "pointer"));
      MAP.on("mouseleave", "zones-fill", () => (MAP.getCanvas().style.cursor = ""));
      zoneLabels(gj);
    });
    const tpp = S.points.tpp3_accident;
    if (tpp) {
      const el = document.createElement("div"); el.className = "tpp-marker"; el.title = "TPP-3 (tank collapse, 29 May 2020)";
      new maplibregl.Marker({ element: el }).setLngLat([tpp[1], tpp[0]])
        .setPopup(new maplibregl.Popup({ offset: 12 }).setHTML("<b>TPP-3</b><br>Tank no. 5 collapsed on 29 May 2020")).addTo(MAP);
    }
    if (ordered.length) setLayer(ordered[0]);
  });

  seg.addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) setLayer(b.dataset.k); });
  document.getElementById("base-seg").addEventListener("click", (e) => {
    const b = e.target.closest("button"); if (!b) return;
    document.querySelectorAll("#base-seg button").forEach((x) => x.setAttribute("aria-pressed", x === b));
    MAP.setLayoutProperty("sat", "visibility", b.dataset.base === "sat" ? "visible" : "none");
    MAP.setLayoutProperty("dark", "visibility", b.dataset.base === "dark" ? "visible" : "none");
  });
  document.getElementById("opacity").addEventListener("change", (e) => {
    if (CURRENT) MAP.setLayoutProperty(CURRENT, "visibility", e.target.checked ? "visible" : "none");
  });
  const leg = document.getElementById("zone-legend");
  leg.innerHTML = ZONES.map((z) => `<li data-z="${z}"><span class="swz" style="color:${ZCOL[z]}"></span>${ZLABEL[z]}</li>`).join("");
  leg.addEventListener("click", (e) => { const li = e.target.closest("li"); if (li) selectZone(li.dataset.z, true); });
}

function zoneLabels(gj) {
  gj.features.forEach((f) => {
    const z = f.properties.zone;
    if (z === "urban") return;
    let best = null, bestA = 0;
    const polys = f.geometry.type === "Polygon" ? [f.geometry.coordinates] : f.geometry.coordinates;
    polys.forEach((p) => { const a = Math.abs(d3.polygonArea(p[0])); if (a > bestA) { bestA = a; best = p[0]; } });
    if (!best) return;
    const c = d3.polygonCentroid(best);
    const el = document.createElement("div"); el.className = "map-label"; el.textContent = ZLABEL[z];
    new maplibregl.Marker({ element: el, anchor: "center" }).setLngLat(c).addTo(MAP);
  });
}

function setLayer(k) {
  const S = window.SITE;
  Object.keys(S.overlays).forEach((id) => MAP.getLayer(id) && MAP.setLayoutProperty(id, "visibility", id === k && document.getElementById("opacity").checked ? "visible" : "none"));
  CURRENT = k;
  document.querySelectorAll("#layer-seg button").forEach((b) => b.setAttribute("aria-checked", b.dataset.k === k));
  const o = S.overlays[k];
  document.getElementById("layer-detail").textContent = `${o.label}: ${o.detail}`;
  const east = k.endsWith("east");
  document.getElementById("cbar-ticks").innerHTML = east
    ? `<span>−${o.lim_mm} mm west</span><span>0</span><span>+${o.lim_mm} mm east</span>`
    : `<span>−${o.lim_mm} mm (sinks)</span><span>0</span><span>+${o.lim_mm} mm (rises)</span>`;
}

function selectZone(z, fly) {
  const S = window.SITE;
  if (MAP.getLayer("zones-hl")) MAP.setFilter("zones-hl", ["==", "zone", z]);
  if (fly) {
    const src = MAP.getSource("zones");
    src && src.getData().then?.((gj) => {
      const f = gj.features.find((ff) => ff.properties.zone === z);
      if (f) { const b = d3.geoBounds(f); MAP.fitBounds(b, { padding: 60, duration: 900 }); }
    });
  }
  const card = document.getElementById("zone-card");
  card.hidden = false;
  const rows = S.timeseries.filter((r) => r.zone === z);
  const fits = S.fits.filter((f) => f.zone === z);
  card.innerHTML = `<h4><span class="swz" style="color:${ZCOL[z]};width:12px;height:12px;border:2px solid;border-radius:3px;display:inline-block"></span>${ZLABEL[z]}</h4>` +
    `<div class="zc-sub">${fits.length ? fits.map((f) => `${f.summer}: ${fmt(f.seasonal_disp_mm)} mm over the summer`).join(" · ") : "Sentinel-1 time series pending"}</div>` +
    `<div id="zc-chart"></div>`;
  if (rows.length) miniSeries(document.getElementById("zc-chart"), rows, z);
}

function miniSeries(el, rows, z) {
  const W = 280, H = 150, m = { l: 34, r: 8, t: 8, b: 22 };
  const byYear = d3.group(rows, (r) => r.summer);
  const doy = (d) => { const t = new Date(d); return (t - new Date(t.getFullYear(), 0, 0)) / 864e5; };
  const x = d3.scaleLinear().domain(d3.extent(rows, (r) => doy(r.date))).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain(d3.extent(rows.map((r) => r.disp_mm).concat([0]))).nice().range([H - m.b, m.t]);
  const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("class", "chart");
  svg.append("g").attr("class", "grid").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4).tickSize(-(W - m.l - m.r)));
  svg.append("line").attr("class", "zero").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(0)).attr("y2", y(0));
  const months = [[152, "Jun"], [182, "Jul"], [213, "Aug"], [244, "Sep"], [274, "Oct"]];
  months.forEach(([d, l]) => { if (d >= x.domain()[0] - 5 && d <= x.domain()[1] + 5) svg.append("text").attr("x", x(d)).attr("y", H - 6).attr("class", "ax-title").attr("text-anchor", "middle").text(l); });
  const dash = { 2020: null, 2025: "5 4" };
  byYear.forEach((rs, yr) => {
    svg.append("path").datum(rs).attr("fill", "none").attr("stroke", ZCOL[z]).attr("stroke-width", 2)
      .attr("stroke-dasharray", dash[yr] || null).attr("d", d3.line().x((r) => x(doy(r.date))).y((r) => y(r.disp_mm)));
    const last = rs[rs.length - 1];
    svg.append("text").attr("x", x(doy(last.date)) - 2).attr("y", y(last.disp_mm) - 6).attr("text-anchor", "end").attr("class", "dlabel").text(yr);
  });
}

/* ---------- results ---------- */
let YEAR, MODE = "rel";
const HIDDEN = new Set();
function results(S) {
  const yrs = [...new Set(S.timeseries.map((r) => r.summer))].sort();
  const tabs = document.getElementById("summer-tabs");
  if (!yrs.length) {
    document.getElementById("chart-scatter").innerHTML = `<p class="note">Sentinel-1 time series are being processed. This chart fills in automatically.</p>`;
    document.getElementById("chart-slopes").innerHTML = `<p class="note">Pending.</p>`;
    return;
  }
  tabs.innerHTML = yrs.map((y, i) => `<button role="tab" aria-selected="${i === 0}" data-y="${y}">Summer ${y}</button>`).join("") +
    `<span class="tab-gap"></span><button class="mode" data-m="rel" aria-pressed="true">Relative to tundra</button><button class="mode" data-m="abs" aria-pressed="false">Absolute vs plateau</button>`;
  tabs.addEventListener("click", (e) => {
    const b = e.target.closest("button"); if (!b) return;
    if (b.dataset.m) { MODE = b.dataset.m; tabs.querySelectorAll(".mode").forEach((x) => x.setAttribute("aria-pressed", x === b)); }
    else { YEAR = +b.dataset.y; tabs.querySelectorAll("button[data-y]").forEach((x) => x.setAttribute("aria-selected", x === b)); }
    scatter(S);
  });
  YEAR = yrs[0];
  scatter(S); slopes(S); fitTable(S); resultText(S);
  document.getElementById("show-table").addEventListener("click", (e) => {
    const t = document.getElementById("fit-table"); t.hidden = !t.hidden; e.target.setAttribute("aria-expanded", !t.hidden);
    e.target.textContent = t.hidden ? "Show the numbers" : "Hide the numbers";
  });
}

function scatter(S) {
  const el = document.getElementById("chart-scatter"); el.innerHTML = "";
  const rel = MODE === "rel";
  const val = (r) => (rel ? r.disp_rel_mm : r.disp_mm);
  const zones = rel ? ZONES.filter((z) => z !== "reference") : ZONES;
  const rows = S.timeseries.filter((r) => r.summer === YEAR && r.tdd_cum != null && val(r) != null);
  const W = 640, H = 400, m = { l: 52, r: 96, t: 14, b: 48 };
  const sq = (r) => Math.sqrt(r.tdd_cum);
  const x = d3.scaleLinear().domain(d3.extent(rows, sq)).nice().range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain(d3.extent(rows.map(val).concat([0]))).nice().range([H - m.b, m.t]);
  const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
  svg.append("g").attr("class", "grid").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSize(-(W - m.l - m.r)));
  svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(7));
  svg.append("line").attr("class", "zero").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(0)).attr("y2", y(0));
  svg.append("text").attr("class", "ax-title").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("√ cumulative thawing degree-days  (√°C·day)");
  svg.append("text").attr("class", "ax-title").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle")
    .text(rel ? "vertical motion minus tundra mean (mm)" : "vertical motion vs plateau (mm, + up)");
  document.getElementById("ch1-sub").textContent = `· summer ${YEAR} · ${rel ? "relative to tundra controls" : "absolute, atmosphere-dominated"}`;
  const labels = [];
  zones.forEach((z) => {
    if (HIDDEN.has(z)) return;
    const rs = rows.filter((r) => r.zone === z).sort((a, b) => sq(a) - sq(b));
    if (!rs.length) return;
    const g = svg.append("g");
    g.append("path").datum(rs).attr("fill", "none").attr("stroke", ZCOL[z]).attr("stroke-width", 2).attr("opacity", 0.85)
      .attr("d", d3.line().x((r) => x(sq(r))).y((r) => y(val(r))));
    const last = rs[rs.length - 1];
    labels.push({ z, y: y(val(last)) });
    g.selectAll("circle").data(rs).join("circle").attr("cx", (r) => x(sq(r))).attr("cy", (r) => y(val(r))).attr("r", 4.5)
      .attr("fill", ZCOL[z]).attr("stroke", css("--card")).attr("stroke-width", 2)
      .on("mousemove", (ev, r) => showTip(`<b>${ZLABEL[z]}</b> · ${r.date}<div class="tt-row"><span class="tt-k">${rel ? "minus tundra" : "vs plateau"}</span><span>${fmt(val(r))} mm</span></div><div class="tt-row"><span class="tt-k">absolute</span><span>${fmt(r.disp_mm)} mm</span></div><div class="tt-row"><span class="tt-k">TDD</span><span>${fmt(r.tdd_cum, 0)} °C·day</span></div><div class="tt-row"><span class="tt-k">pixels</span><span>${r.n_px.toLocaleString()}</span></div>`, ev))
      .on("mouseleave", hideTip);
  });
  labels.sort((a, b) => a.y - b.y);
  for (let i = 1; i < labels.length; i++) if (labels[i].y - labels[i - 1].y < 14) labels[i].y = labels[i - 1].y + 14;
  labels.forEach((l) => svg.append("text").attr("class", "dlabel").attr("x", W - m.r + 8).attr("y", l.y + 4).text(ZLABEL[l.z]));
  const lg = document.getElementById("ch1-legend");
  lg.innerHTML = zones.map((z) => `<span data-z="${z}" class="${HIDDEN.has(z) ? "off" : ""}"><i style="background:${ZCOL[z]}"></i>${ZLABEL[z]}</span>`).join("");
  lg.onclick = (e) => { const s2 = e.target.closest("span"); if (!s2) return; const z = s2.dataset.z; HIDDEN.has(z) ? HIDDEN.delete(z) : HIDDEN.add(z); scatter(S); };
}

function slopes(S) {
  const el = document.getElementById("chart-slopes"); el.innerHTML = "";
  const want = [["industrial", "controls", "Industrial − tundra"], ["accident", "controls", "TPP-3 − tundra"],
    ["urban", "controls", "City − tundra"], ["control_NW", "control_SE", "Tundra NW − SE"]];
  const yrs = [...new Set(S.comparisons.map((c) => c.summer))].sort();
  const rowsC = want.map(([a, b, l]) => ({ a, b, l, by: Object.fromEntries(yrs.map((yr) => [yr, S.comparisons.find((c) => c.zone_a === a && c.zone_b === b && c.summer === yr)])) }));
  const W = 460, rowH = 54, m = { l: 132, r: 24, t: 12, b: 70 };
  const H = m.t + m.b + rowH * rowsC.length;
  const vals = S.comparisons.flatMap((c) => [c.diff_b - 2 * c.se, c.diff_b + 2 * c.se]).filter(isFinite).concat([0]);
  const x = d3.scaleLinear().domain(d3.extent(vals)).nice().range([m.l, W - m.r]);
  const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
  svg.append("g").attr("class", "grid").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5).tickSize(-(H - m.t - m.b)));
  svg.append("line").attr("class", "zero").attr("x1", x(0)).attr("x2", x(0)).attr("y1", m.t).attr("y2", H - m.b);
  svg.append("text").attr("class", "ax-title").attr("x", (m.l + W - m.r) / 2).attr("y", H - m.b + 38).attr("text-anchor", "middle").text("Δ slope, mm per √(°C·day)   ← first sinks more");
  const off = { [yrs[0]]: -9, [yrs[1]]: 9 };
  rowsC.forEach((rw, i) => {
    const cy = m.t + rowH * i + rowH / 2;
    svg.append("text").attr("x", m.l - 12).attr("y", cy + 4).attr("text-anchor", "end").attr("fill", css("--ink-2")).attr("font-size", 13).text(rw.l);
    const col = rw.a === "control_NW" ? ZCOL.control_SE : ZCOL[rw.a];
    yrs.forEach((yr) => {
      const c = rw.by[yr]; if (!c || c.diff_b == null) return;
      const yy = cy + (off[yr] || 0);
      svg.append("line").attr("x1", x(c.diff_b - 2 * c.se)).attr("x2", x(c.diff_b + 2 * c.se)).attr("y1", yy).attr("y2", yy).attr("stroke", col).attr("stroke-width", 2).attr("opacity", 0.55);
      const first = yr === yrs[0];
      const mk = first ? svg.append("circle").attr("r", 5).attr("cx", x(c.diff_b)).attr("cy", yy)
        : svg.append("rect").attr("width", 9).attr("height", 9).attr("rx", 2).attr("x", x(c.diff_b) - 4.5).attr("y", yy - 4.5);
      mk.attr("fill", col).attr("stroke", css("--card")).attr("stroke-width", 2)
        .on("mousemove", (ev) => showTip(`<b>${rw.l}</b> · ${yr}<div class="tt-row"><span class="tt-k">Δ slope</span><span>${fmt(c.diff_b, 2)} ± ${fmt(2 * c.se, 2)}</span></div><div class="tt-row"><span class="tt-k">z</span><span>${fmt(c.z, 1)}</span></div><div class="tt-row"><span class="tt-k">difference at summer end</span><span>${fmt(c.seasonal_diff_mm)} mm</span></div>`, ev))
        .on("mouseleave", hideTip);
    });
  });
  const lg = svg.append("g").attr("transform", `translate(${m.l},${H - 8})`);
  yrs.forEach((yr, i) => {
    const g = lg.append("g").attr("transform", `translate(${i * 80},-4)`);
    i === 0 ? g.append("circle").attr("r", 4.5).attr("fill", css("--ink-2")) : g.append("rect").attr("x", -4).attr("y", -4).attr("width", 8).attr("height", 8).attr("rx", 2).attr("fill", css("--ink-2"));
    g.append("text").attr("x", 10).attr("y", 4).attr("fill", css("--muted")).attr("font-size", 12).text(yr);
  });
}

function fitTable(S) {
  const rows = S.fits.map((f) => `<tr><td>${ZLABEL[f.zone] || f.zone}</td><td>${f.summer}</td><td>${fmt(f.b, 2)}</td><td>${fmt(f.se_b, 2)}</td><td>${fmt(f.r2, 2)}</td><td>${f.n}</td><td>${fmt(f.seasonal_disp_mm)}</td></tr>`).join("");
  const name = (z) => (z === "controls" ? "tundra mean" : ZLABEL[z]);
  const crow = S.comparisons.map((c) => `<tr><td>${name(c.zone_a)} − ${name(c.zone_b)}</td><td>${c.summer}</td><td>${fmt(c.diff_b, 2)}</td><td>${fmt(c.se, 2)}</td><td>${fmt(c.z, 1)}</td><td>${fmt(c.seasonal_diff_mm)}</td></tr>`).join("");
  document.getElementById("fit-table").innerHTML =
    `<p class="card-note">Paired comparisons (robust to shared atmosphere)</p><table><thead><tr><th>Comparison</th><th>Summer</th><th>Δb</th><th>σ</th><th>z</th><th>end-of-summer Δ (mm)</th></tr></thead><tbody>${crow}</tbody></table>` +
    `<p class="card-note" style="margin-top:18px">Absolute fits vs plateau (atmosphere-dominated)</p><table><thead><tr><th>Zone</th><th>Summer</th><th>b (mm/√°C·d)</th><th>σ</th><th>R²</th><th>dates</th><th>summer total (mm)</th></tr></thead><tbody>${rows}</tbody></table>`;
}

function resultText(S) {
  const V = S.verdict; if (!V || V.status !== "ready") return;
  const parts = Object.entries(V.summers).map(([yr, v]) => {
    const d = v.seasonal_industrial_minus_controls_mm;
    return `${yr}: industrial ground ended the summer ${fmt(Math.abs(d))} mm ${d < 0 ? "lower" : "higher"} than the tundra (z = ${fmt(v.z_industrial_vs_controls, 1)}).`;
  });
  document.getElementById("result-title").textContent = "Industry and tundra move together, and so does the air";
  document.getElementById("result-lede").textContent =
    "All lowland zones rise and fall together by tens of millimetres from one pass to the next: that common motion is water vapour between the lowland and the plateau reference, not ground. " +
    "Subtracting the tundra removes it and leaves the question that matters. " + parts.join(" ");
}

/* ---------- local, reproducible signal ---------- */
const LC_ORDER = ["herbaceous wetland", "moss and lichen", "grassland", "water edge", "tree cover", "bare / sparse", "built-up"];
function local(S) {
  const L = S.local; if (!L) { document.getElementById("local").hidden = true; return; }
  const wet = L.landcover.find((r) => r.code === 90);
  document.getElementById("local-stats").innerHTML = [
    { v: fmt(L.r_local, 2), k: `pixel correlation of local thaw motion, 2020 vs 2025 (${L.n_independent.toLocaleString()} samples 150 m apart)` },
    { v: `${fmt(L.r_null_shift_mean, 2)}<small> ± ${fmt(L.r_null_shift_std, 2)}</small>`, k: "same correlation when the 2025 map is shifted 6 km: what chance looks like" },
    { v: `${fmt(L.hotspot_px_both / L.hotspot_px_expected_if_independent, 1)}×`, k: `more pixels sink by > 10 mm in both summers than if the summers were unrelated (${fmt(L.hotspot_area_km2, 0)} km²)` },
    { v: `${fmt(wet.local_2020_mm)} / ${fmt(wet.local_2025_mm)}<small>mm</small>`, k: "wetland ground vs its 2 km surroundings, end of summer 2020 / 2025", hl: true }
  ].map((t) => `<div class="stat ${t.hl ? "hl" : ""}"><div class="v">${t.v}</div><div class="k">${t.k}</div></div>`).join("");
  localSeries(S); localLandcover(S); localFacilities(S);
}

function localSeries(S) {
  const L = S.local, el = document.getElementById("chart-local-series"); el.innerHTML = "";
  const yrs = Object.keys(L.series);
  const W = 640, H = 360, m = { l: 48, r: 132, t: 30, b: 46 }, gap = 26;
  const pw = (W - m.l - m.r - gap) / yrs.length;
  const all = yrs.flatMap((y) => Object.values(L.series[y].by_class).flatMap((v) => v.map((d) => d - v[0])));
  const y = d3.scaleLinear().domain(d3.extent(all.concat([0]))).nice().range([H - m.b, m.t]);
  const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
  svg.append("text").attr("class", "ax-title").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("motion vs 2 km surroundings (mm, + up)");
  const endLabels = [];
  yrs.forEach((yr, i) => {
    const s = L.series[yr], x0 = m.l + i * (pw + gap);
    const x = d3.scaleLinear().domain(d3.extent(s.sqrt_tdd)).nice().range([x0, x0 + pw]);
    const g = svg.append("g");
    g.append("g").attr("class", "grid").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).ticks(5).tickSize(-pw).tickFormat(i ? "" : null));
    g.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(4));
    g.append("line").attr("class", "zero").attr("x1", x0).attr("x2", x0 + pw).attr("y1", y(0)).attr("y2", y(0));
    g.append("text").attr("x", x0).attr("y", 16).attr("fill", css("--ink")).attr("font-family", "Space Grotesk").attr("font-weight", 600).attr("font-size", 14).text(`Summer ${yr}`);
    const cls = LC_ORDER.filter((c) => s.by_class[c]);
    cls.slice().reverse().forEach((c) => {
      const v = s.by_class[c].map((d) => d - s.by_class[c][0]), hl = c === "herbaceous wetland";
      const pts = v.map((d, k) => [s.sqrt_tdd[k], d, s.dates[k]]);
      g.append("path").datum(pts).attr("fill", "none").attr("stroke", hl ? css("--sub") : "#5b5f66").attr("stroke-width", hl ? 2.5 : 1.5).attr("opacity", hl ? 1 : 0.8)
        .attr("d", d3.line().x((p) => x(p[0])).y((p) => y(p[1])));
      g.selectAll(null).data(pts).join("circle").attr("cx", (p) => x(p[0])).attr("cy", (p) => y(p[1])).attr("r", hl ? 4.5 : 3)
        .attr("fill", hl ? css("--sub") : "#5b5f66").attr("stroke", css("--card")).attr("stroke-width", 2)
        .on("mousemove", (ev, p) => showTip(`<b>${c}</b> · ${p[2].slice(0, 4)}-${p[2].slice(4, 6)}-${p[2].slice(6)}<div class="tt-row"><span class="tt-k">vs surroundings</span><span>${fmt(p[1])} mm</span></div><div class="tt-row"><span class="tt-k">√TDD</span><span>${fmt(p[0])}</span></div>`, ev))
        .on("mouseleave", hideTip);
      if (i === yrs.length - 1) endLabels.push({ c, hl, y: y(v[v.length - 1]) });
    });
  });
  endLabels.sort((a, b) => a.y - b.y);
  for (let k = 1; k < endLabels.length; k++) if (endLabels[k].y - endLabels[k - 1].y < 13) endLabels[k].y = endLabels[k - 1].y + 13;
  endLabels.forEach((l) => svg.append("text").attr("class", "dlabel").attr("x", W - m.r + 8).attr("y", l.y + 4)
    .attr("fill", l.hl ? css("--ink") : css("--muted")).attr("font-weight", l.hl ? 600 : 400).text(l.c));
  svg.append("text").attr("class", "ax-title").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("√ cumulative thawing degree-days  (√°C·day)");
}

function localLandcover(S) {
  const L = S.local, el = document.getElementById("chart-local-lc"); el.innerHTML = "";
  const rows = LC_ORDER.map((c) => L.landcover.find((r) => r.class === c)).filter(Boolean);
  const W = 460, rowH = 40, m = { l: 140, r: 20, t: 10, b: 62 }, H = m.t + m.b + rowH * rows.length;
  const vals = rows.flatMap((r) => [r.local_2020_mm - 2 * r.se_2020, r.local_2020_mm + 2 * r.se_2020, r.local_2025_mm - 2 * r.se_2025, r.local_2025_mm + 2 * r.se_2025]).concat([0]);
  const x = d3.scaleLinear().domain(d3.extent(vals)).nice().range([m.l, W - m.r]);
  const svg = d3.select(el).append("svg").attr("viewBox", `0 0 ${W} ${H}`);
  svg.append("g").attr("class", "grid").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5).tickSize(-(H - m.t - m.b)));
  svg.append("line").attr("class", "zero").attr("x1", x(0)).attr("x2", x(0)).attr("y1", m.t).attr("y2", H - m.b);
  svg.append("text").attr("class", "ax-title").attr("x", (m.l + W - m.r) / 2).attr("y", H - m.b + 36).attr("text-anchor", "middle").text("end-of-summer motion vs surroundings (mm)  ← sinks more");
  rows.forEach((r, i) => {
    const cy = m.t + rowH * i + rowH / 2, hl = r.code === 90, col = hl ? css("--sub") : css("--ink-2");
    svg.append("text").attr("x", m.l - 12).attr("y", cy + 4).attr("text-anchor", "end").attr("fill", hl ? css("--ink") : css("--ink-2")).attr("font-weight", hl ? 600 : 400).attr("font-size", 13).text(r.class);
    [["2020", -7], ["2025", 7]].forEach(([yr, off], j) => {
      const v = r[`local_${yr}_mm`], se = r[`se_${yr}`], yy = cy + off;
      svg.append("line").attr("x1", x(v - 2 * se)).attr("x2", x(v + 2 * se)).attr("y1", yy).attr("y2", yy).attr("stroke", col).attr("stroke-width", 2).attr("opacity", 0.5);
      const mk = j === 0 ? svg.append("circle").attr("r", 5).attr("cx", x(v)).attr("cy", yy)
        : svg.append("rect").attr("width", 9).attr("height", 9).attr("rx", 2).attr("x", x(v) - 4.5).attr("y", yy - 4.5);
      mk.attr("fill", col).attr("stroke", css("--card")).attr("stroke-width", 2)
        .on("mousemove", (ev) => showTip(`<b>${r.class}</b> · ${yr}<div class="tt-row"><span class="tt-k">vs surroundings</span><span>${fmt(v)} ± ${fmt(2 * se)} mm</span></div><div class="tt-row"><span class="tt-k">area</span><span>${fmt(r.area_km2, 0)} km²</span></div><div class="tt-row"><span class="tt-k">sinks > 10 mm both summers</span><span>${fmt(100 * r.frac_hotspot)} %</span></div>`, ev))
        .on("mouseleave", hideTip);
    });
  });
  const lg = svg.append("g").attr("transform", `translate(${m.l},${H - 8})`);
  lg.append("circle").attr("cy", -4).attr("r", 4.5).attr("fill", css("--ink-2"));
  lg.append("text").attr("x", 10).attr("y", 0).attr("fill", css("--muted")).attr("font-size", 12).text("2020");
  lg.append("rect").attr("x", 76).attr("y", -8).attr("width", 8).attr("height", 8).attr("rx", 2).attr("fill", css("--ink-2"));
  lg.append("text").attr("x", 90).attr("y", 0).attr("fill", css("--muted")).attr("font-size", 12).text("2025   (±2σ)");
}

function localFacilities(S) {
  const L = S.local, el = document.getElementById("local-fac");
  const rows = L.facilities.slice().sort((a, b) => a.d_mean - b.d_mean);
  const lu = Object.fromEntries(L.facilities_by_landuse.map((r) => [r.landuse, r]));
  el.innerHTML = `<table><thead><tr><th>Site</th><th>Type</th><th>2020</th><th>2025</th><th></th></tr></thead><tbody>` +
    rows.map((r) => `<tr data-lon="${r.lon}" data-lat="${r.lat}"><td>${r.label}</td><td class="muted">${r.landuse}</td><td>${r.d_2020 > 0 ? "+" : ""}${fmt(r.d_2020)}</td><td>${r.d_2025 > 0 ? "+" : ""}${fmt(r.d_2025)}</td>` +
      `<td><span class="badge ${r.d_mean < 0 ? "sink" : "rise"}">${r.d_mean < 0 ? "sinks more" : "sinks less"}</span></td></tr>`).join("") + `</tbody></table>` +
    `<p class="card-note">mm, end of summer, inside the site minus a 0.3–1.5 km ring around it. Listed: ${rows.length} of ${L.facilities_all_n} OSM sites where the difference has the same sign and exceeds 2σ in both summers (site-level correlation between summers r = ${fmt(L.facility_r, 2)}). ` +
    `On average, industrial sites: ${fmt(lu.industrial?.d_2020)} / ${fmt(lu.industrial?.d_2025)} mm; landfills and tailings: ${fmt(lu.landfill?.d_2020)} / ${fmt(lu.landfill?.d_2025)} mm. Click a row to find it on the map.</p>`;
  el.querySelectorAll("tr[data-lon]").forEach((tr) => tr.addEventListener("click", () => {
    if (!MAP) return;
    setLayer("local_persistent");
    document.getElementById("explore").scrollIntoView({ behavior: "smooth" });
    MAP.flyTo({ center: [+tr.dataset.lon, +tr.dataset.lat], zoom: 13, duration: 1400 });
  }));
}

/* ---------- verdict ---------- */
function verdict(S) {
  const V = S.verdict, grid = document.getElementById("verdict-grid");
  const H = [
    { k: "climatic", t: "A · Climatic thaw", d: "Both tundra controls subside significantly with √TDD, relative to bedrock." },
    { k: "anthropogenic", t: "B · Anthropogenic thaw", d: "Industrial ground subsides significantly more than the tundra (paired test)." },
    { k: "structural", t: "C · Stable ground", d: "No significant seasonal subsidence in industry or tundra, relative to bedrock." }
  ];
  if (!V || V.status !== "ready") {
    grid.innerHTML = H.map((h) => `<div class="vcard"><div class="vh"><h4>${h.t}</h4><span class="badge no">pending</span></div><div class="meter"><i></i></div><p>${h.d}</p></div>`).join("");
    return;
  }
  const n = V.n_summers;
  const state = (c, k) => c.yes === n ? ["good", "✓ supported", css("--good"), 100]
    : (k === "climatic" && V.climatic_local) ? ["warn", "◐ local evidence", css("--warning"), 50]
    : c.yes > 0 ? ["warn", "◐ partly", css("--warning"), 50]
    : c.untestable === n ? ["no", "? not testable yet", "#3a3f48", 0]
    : ["no", "✕ not supported", "#3a3f48", 0];
  grid.innerHTML = H.map((h) => {
    const c = V[h.k], [cls, lab, col, w] = state(c, h.k);
    const why = (h.k === "climatic" && V.climatic_local) ? "Absolute motion vs bedrock is still atmosphere-dominated, but ice-rich wetland ground sinks <strong>relative to its surroundings</strong> along the thaw curve in both summers (see “The pattern that repeats”)."
      : c.untestable === n ? "Absolute motion is dominated by atmosphere (σ too large); needs tropospheric correction."
      : `Met in <strong>${c.yes} of ${n}</strong> summers.`;
    return `<div class="vcard"><div class="vh"><h4>${h.t}</h4><span class="badge ${cls}">${lab}</span></div>
      <div class="meter"><i data-w="${Math.max(w, 4)}" style="background:${col}"></i></div><p>${h.d} ${why}</p></div>`;
  }).join("");
  setTimeout(() => grid.querySelectorAll(".meter i").forEach((i) => (i.style.width = i.dataset.w + "%")), 300);
  document.querySelectorAll(".hyp").forEach((card) => {
    const c = V[card.dataset.h], [cls, lab] = state(c, card.dataset.h), el = card.querySelector(".hyp-verdict");
    card.classList.toggle("supported", cls === "good"); card.classList.toggle("partial", cls === "warn");
    el.innerHTML = `<span class="badge ${cls}">${lab}</span>`;
  });
  const S2 = Object.entries(V.summers);
  const zs = S2.map(([yr, v]) => `${yr}: z = ${fmt(v.z_industrial_vs_controls, 1)}, ${fmt(v.seasonal_industrial_minus_controls_mm)} mm`).join("; ");
  const ctl = S2.map(([yr, v]) => `${yr}: ${fmt(Math.abs(v.seasonal_controlNW_minus_SE_mm))} mm`).join(", ");
  let head = V.anthropogenic.yes === n ? "Industrial ground sinks faster than natural tundra."
    : V.anthropogenic.yes === 0 ? "No sign that the city is sinking faster than the tundra around it."
    : "Industrial ground sank faster than tundra in one summer but not the other.";
  document.getElementById("verdict-title").textContent = head;
  document.getElementById("verdict-text").innerHTML =
    `<p>${head} The difference between industrial and natural ground changes sign between summers (${zs}) and stays within 2σ. The two tundra controls differ from each other by ${ctl}: natural variability is as large as any industrial excess.</p>` +
    `<p>What this rules out: a strong, city-wide <em>anthropogenic</em> thaw signal larger than ~1–1.5 cm per summer (2σ) relative to the surrounding tundra. What it cannot yet say: whether all of the lowland, tundra included, is sinking relative to bedrock. That needs the atmosphere removed first.</p>` +
    `<p>Rule: "supported" when a slope differs from zero by more than 2σ. Paired tests subtract the tundra on each date, so shared atmosphere cancels. Absolute tests are marked not testable when σ exceeds ${0.25} mm/√(°C·day).</p>`;
}

/* ---------- NISAR ---------- */
function nisar(S) {
  const N = S.nisar, box = document.getElementById("nisar-stats"), card = document.getElementById("nisar-card");
  if (!N || !N.late) { box.innerHTML = `<p class="note">NISAR results pending.</p>`; return; }
  const w = ["early", "late"].filter((k) => N[k]);
  box.innerHTML = w.map((k) => `<div class="stat"><div class="v">${fmt(N[k].flat_abs_east_p50_mm)}<small>mm</small></div><div class="k">median |east-west| motion on flat ground, ${k} summer (12 days)</div></div>
    <div class="stat"><div class="v">${fmt(N[k].flat_abs_up_p50_mm)}<small>mm</small></div><div class="k">median |vertical| motion on flat ground, ${k} summer (12 days)</div></div>`).join("");
  // per-zone grouped dot plot: up vs east, late window
  const k = "late", pz = N[k].per_zone, zs = ZONES.filter((z) => pz[z]);
  const W = 460, rowH = 40, m = { l: 112, r: 20, t: 34, b: 40 }, H = m.t + m.b + rowH * zs.length;
  const vals = zs.flatMap((z) => [pz[z].up_mean_mm, pz[z].east_mean_mm]).concat([0]);
  const x = d3.scaleLinear().domain(d3.extent(vals)).nice().range([m.l, W - m.r]);
  const svg = d3.select(card).append("svg").attr("viewBox", `0 0 ${W} ${H}`).attr("class", "chart");
  svg.append("text").attr("x", 0).attr("y", 14).attr("fill", css("--ink")).attr("font-family", "Space Grotesk").attr("font-weight", 600).attr("font-size", 15)
    .text("Zone means · late summer 2026 (12-day pairs)");
  svg.append("g").attr("class", "grid").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5).tickSize(-(H - m.t - m.b)));
  svg.append("line").attr("class", "zero").attr("x1", x(0)).attr("x2", x(0)).attr("y1", m.t).attr("y2", H - m.b);
  zs.forEach((z, i) => {
    const cy = m.t + rowH * i + rowH / 2;
    svg.append("text").attr("x", m.l - 12).attr("y", cy + 4).attr("text-anchor", "end").attr("fill", css("--ink-2")).attr("font-size", 13).text(ZLABEL[z]);
    svg.append("line").attr("x1", x(pz[z].up_mean_mm)).attr("x2", x(pz[z].east_mean_mm)).attr("y1", cy).attr("y2", cy).attr("stroke", "#3a3f48").attr("stroke-width", 2);
    svg.append("circle").attr("cx", x(pz[z].up_mean_mm)).attr("cy", cy).attr("r", 5.5).attr("fill", ZCOL[z]).attr("stroke", css("--card")).attr("stroke-width", 2)
      .on("mousemove", (ev) => showTip(`<b>${ZLABEL[z]}</b><div class="tt-row"><span class="tt-k">vertical</span><span>${fmt(pz[z].up_mean_mm)} mm</span></div><div class="tt-row"><span class="tt-k">east</span><span>${fmt(pz[z].east_mean_mm)} mm</span></div>`, ev)).on("mouseleave", hideTip);
    svg.append("rect").attr("x", x(pz[z].east_mean_mm) - 5).attr("y", cy - 5).attr("width", 10).attr("height", 10).attr("rx", 2)
      .attr("fill", "none").attr("stroke", ZCOL[z]).attr("stroke-width", 2);
  });
  const lg = svg.append("g").attr("transform", `translate(${m.l},${H - 8})`);
  lg.append("circle").attr("cx", 0).attr("cy", -4).attr("r", 4.5).attr("fill", css("--ink-2"));
  lg.append("text").attr("x", 10).attr("y", 0).attr("fill", css("--muted")).attr("font-size", 12).text("vertical");
  lg.append("rect").attr("x", 80).attr("y", -9).attr("width", 9).attr("height", 9).attr("rx", 2).attr("fill", "none").attr("stroke", css("--ink-2")).attr("stroke-width", 2);
  lg.append("text").attr("x", 96).attr("y", 0).attr("fill", css("--muted")).attr("font-size", 12).text("east-west   (mm, 12 days)");
  const e = N.late.flat_abs_east_p50_mm, u = N.late.flat_abs_up_p50_mm;
  const p = document.createElement("p"); p.className = "card-note";
  p.textContent = e >= 0.6 * u
    ? "On flat ground the east-west component is as large as the vertical one. With single 12-day pairs this mostly reflects residual atmospheric noise, so NISAR cannot yet confirm that the motion is purely vertical: stacking more pairs is the next step."
    : "On flat ground the east-west component is small compared with the vertical one, supporting the vertical assumption used for Sentinel-1.";
  card.appendChild(p);
}

/* ---------- scroll reveal ---------- */
function reveal() {
  const els = document.querySelectorAll(".section .wrap > *, .hyp, .step, .vcard, .limit");
  els.forEach((e) => e.classList.add("reveal"));
  const io = new IntersectionObserver((es) => es.forEach((en) => { if (en.isIntersecting) { en.target.classList.add("in"); io.unobserve(en.target); } }), { threshold: 0.12 });
  els.forEach((e) => io.observe(e));
}
