/* Where ORCHID data files are served from. On orchid-project.org they come from the separate data site
   (repo orchid-calphad/orchid-data); elsewhere (local preview, tests) from the sibling ../data/ folder. */
(function () {
  var h = location.hostname;
  window.ORCHID_DATA = /(^|\.)orchid-project\.org$/.test(h) && h.indexOf("data.") !== 0 ? "https://data.orchid-project.org/" : "../data/";
  window.dataURL = function (p) { return window.ORCHID_DATA + p; };
  document.addEventListener("DOMContentLoaded", function () {
    if (window.ORCHID_DATA === "../data/") return;
    document.querySelectorAll('a[href^="../data/literature/"],a[href^="../data/phase-diagrams/"],a[href^="../data/hulls/"]').forEach(function (a) {
      a.setAttribute("href", window.ORCHID_DATA + a.getAttribute("href").slice("../data/".length));
    });
  });
})();
