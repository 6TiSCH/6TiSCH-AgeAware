(function (t) {
  function e(e) {
    for (
      var a, o, r = e[0], l = e[1], u = e[2], d = 0, h = [];
      d < r.length;
      d++
    )
      (o = r[d]), s[o] && h.push(s[o][0]), (s[o] = 0);
    for (a in l) Object.prototype.hasOwnProperty.call(l, a) && (t[a] = l[a]);
    c && c(e);
    while (h.length) h.shift()();
    return i.push.apply(i, u || []), n();
  }
  function n() {
    for (var t, e = 0; e < i.length; e++) {
      for (var n = i[e], a = !0, r = 1; r < n.length; r++) {
        var l = n[r];
        0 !== s[l] && (a = !1);
      }
      a && (i.splice(e--, 1), (t = o((o.s = n[0]))));
    }
    return t;
  }
  var a = {},
    s = { app: 0 },
    i = [];
  function o(e) {
    if (a[e]) return a[e].exports;
    var n = (a[e] = { i: e, l: !1, exports: {} });
    return t[e].call(n.exports, n, n.exports, o), (n.l = !0), n.exports;
  }
  (o.m = t),
    (o.c = a),
    (o.d = function (t, e, n) {
      o.o(t, e) || Object.defineProperty(t, e, { enumerable: !0, get: n });
    }),
    (o.r = function (t) {
      "undefined" !== typeof Symbol &&
        Symbol.toStringTag &&
        Object.defineProperty(t, Symbol.toStringTag, { value: "Module" }),
        Object.defineProperty(t, "__esModule", { value: !0 });
    }),
    (o.t = function (t, e) {
      if ((1 & e && (t = o(t)), 8 & e)) return t;
      if (4 & e && "object" === typeof t && t && t.__esModule) return t;
      var n = Object.create(null);
      if (
        (o.r(n),
        Object.defineProperty(n, "default", { enumerable: !0, value: t }),
        2 & e && "string" != typeof t)
      )
        for (var a in t)
          o.d(
            n,
            a,
            function (e) {
              return t[e];
            }.bind(null, a)
          );
      return n;
    }),
    (o.n = function (t) {
      var e =
        t && t.__esModule
          ? function () {
              return t["default"];
            }
          : function () {
              return t;
            };
      return o.d(e, "a", e), e;
    }),
    (o.o = function (t, e) {
      return Object.prototype.hasOwnProperty.call(t, e);
    }),
    (o.p = "/");
  var r = (window["webpackJsonp"] = window["webpackJsonp"] || []),
    l = r.push.bind(r);
  (r.push = e), (r = r.slice());
  for (var u = 0; u < r.length; u++) e(r[u]);
  var c = l;
  i.push([0, "chunk-vendors"]), n();
})({
  0: function (t, e, n) {
    t.exports = n("56d7");
  },
  "56d7": function (t, e, n) {
    "use strict";
    n.r(e);
    n("cadf"), n("551c"), n("097d");
    var a = n("2b0e"),
      s = n("bb71");
    n("da64");
    a["a"].use(s["a"], { iconfont: "md" });
    var i = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-app",
          [
            n("TheToolbar"),
            n("TheNavigationDrawer"),
            n(
              "v-content",
              { attrs: { fluid: "" } },
              [
                n(
                  "v-layout",
                  {
                    attrs: {
                      "align-center": "",
                      "justify-center": "",
                      row: "",
                      wrap: "",
                      "fill-height": "",
                    },
                  },
                  [n("router-view")],
                  1
                ),
              ],
              1
            ),
            n("TheFooter"),
            n("TheAboutDialog"),
            n("TheSettingsDialog"),
          ],
          1
        );
      },
      o = [],
      r = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-toolbar",
          { attrs: { dense: "", app: "" } },
          [
            n("v-toolbar-title", { staticClass: "headline text-uppercase" }, [
              n("span", [t._v("6TiSCH")]),
              n("span", { staticClass: "font-weight-light" }, [
                t._v("Simulator"),
              ]),
            ]),
            n("v-spacer"),
            n("TheConnectionStatus"),
            n("v-toolbar-side-icon", {
              attrs: { disabled: t.disabled },
              on: {
                click: function (e) {
                  t.$_app_toggleNavigationDrawer();
                },
              },
            }),
          ],
          1
        );
      },
      l = [],
      u = {
        computed: {
          $_app_status: {
            get: function () {
              return this.$store.getters["app/status"];
            },
            set: function (t) {
              this.$store.dispatch("app/setStatus", t);
            },
          },
          $_app_shutdownDialog: {
            get: function () {
              return this.$store.getters["app/shutdownDialog"];
            },
            set: function (t) {
              this.$store.dispatch("app/setShutdownDialog", t);
            },
          },
          $_app_autoReload: {
            get: function () {
              return this.$store.getters["app/autoReload"];
            },
            set: function (t) {
              this.$store.dispatch("app/setAutoReload", t);
            },
          },
          $_app_aboutDialog: {
            get: function () {
              return this.$store.getters["app/aboutDialog"];
            },
            set: function (t) {
              this.$store.dispatch("app/setAboutDialog", t);
            },
          },
          $_app_settingsDialog: {
            get: function () {
              return this.$store.getters["app/settingsDialog"];
            },
            set: function (t) {
              this.$store.dispatch("app/setConfigDialog", t);
            },
          },
        },
        methods: {
          $_app_toggleNavigationDrawer: function () {
            this.$store.dispatch("app/toggleNavigationDrawer");
          },
        },
      },
      c = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "span",
          [
            n(
              "v-tooltip",
              { attrs: { bottom: "" } },
              [
                n(
                  "v-icon",
                  { attrs: { slot: "activator" }, slot: "activator" },
                  [
                    t._v(
                      "\n      " +
                        t._s(t._f("icon")(t.$_simulator_connectionStatus)) +
                        "\n    "
                    ),
                  ]
                ),
                n("span", [t._v(t._s(t.$_simulator_connectionStatus))]),
              ],
              1
            ),
            n(
              "v-dialog",
              {
                attrs: { persistent: "", "max-width": "400" },
                model: {
                  value: t.dialog,
                  callback: function (e) {
                    t.dialog = e;
                  },
                  expression: "dialog",
                },
              },
              [
                n(
                  "v-card",
                  [
                    n(
                      "v-layout",
                      [
                        n(
                          "v-card-title",
                          {
                            staticClass: "title",
                            attrs: { "primary-title": "" },
                          },
                          [
                            t._v(
                              "\n          Disconnected from Backend\n        "
                            ),
                          ]
                        ),
                      ],
                      1
                    ),
                    n(
                      "v-layout",
                      [
                        n(
                          "v-flex",
                          [
                            n("v-card-text", { staticClass: "subheading" }, [
                              t._v(
                                "\n            Reloading the WebApp in " +
                                  t._s(t.remaining) +
                                  "s ...\n          "
                              ),
                            ]),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                    n(
                      "v-layout",
                      {
                        attrs: {
                          "align-center": "",
                          "justify-center": "",
                          row: "",
                        },
                      },
                      [
                        n(
                          "v-card-actions",
                          [
                            n(
                              "v-btn",
                              {
                                attrs: { color: "error" },
                                on: {
                                  click: function (e) {
                                    e.stopPropagation(), t.stopTimer();
                                  },
                                },
                              },
                              [t._v("\n            Cancel\n          ")]
                            ),
                            n(
                              "v-btn",
                              {
                                on: {
                                  click: function (e) {
                                    e.stopPropagation(), t.reload();
                                  },
                                },
                              },
                              [t._v("\n            Reload Now\n          ")]
                            ),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      d = [],
      h = {
        computed: {
          $_simulator_connectionStatus: function () {
            return this.$store.getters["simulator/connectionStatus"];
          },
          $_simulator_operationalStatus: function () {
            return this.$store.getters["simulator/operationalStatus"];
          },
          $_simulator_defaultSettings: {
            get: function () {
              return this.$store.getters["simulator/defaultSettings"];
            },
            set: function (t) {
              this.$store.dispatch("simulator/saveDefaultSettings", t);
            },
          },
          $_simulator_runningSettings: {
            get: function () {
              return this.$store.getters["simulator/runningSettings"];
            },
            set: function (t) {
              this.$store.dispatch("simulator/updateRunningSettings", t);
            },
          },
          $_simulator_availableSFs: function () {
            return this.$store.getters["simulator/availableSFs"];
          },
          $_simulator_availableConnectivities: function () {
            return this.$store.getters["simulator/availableConnectivities"];
          },
          $_simulator_availableTraceFiles: function () {
            return this.$store.getters["simulator/availableTraceFiles"];
          },
          $_simulator_gitInfo: function () {
            return this.$store.getters["simulator/gitInfo"];
          },
          $_simulator_crashReport: function () {
            return this.$store.getters["simulator/crashReport"];
          },
        },
        methods: {
          $_simulator_start: function () {
            return this.$store.dispatch("simulator/start", this.$_eel);
          },
          $_simulator_pause: function () {
            return this.$store.dispatch("simulator/pause", this.$_eel);
          },
          $_simulator_resume: function () {
            return this.$store.dispatch("simulator/resume", this.$_eel);
          },
          $_simulator_abort: function () {
            return this.$store.dispatch("simulator/abort", this.$_eel);
          },
          $_simulator_clearCrashReport: function () {
            this.$store.dispatch("simulator/clearCrashReport");
          },
        },
      },
      m = {
        filters: {
          icon: function (t) {
            return "connected" === t ? "cloud_done" : "cloud_off";
          },
        },
        mixins: [u, h],
        data: function () {
          return { dialog: !1, timeToReload: 10, timePassed: 0 };
        },
        computed: {
          remaining: function () {
            return this.timeToReload - this.timePassed;
          },
        },
        watch: {
          $_simulator_connectionStatus: function (t) {
            if ("disconnected" === t && !0 === this.$_app_autoReload) {
              (this.$_app_settingsDialog = !1),
                (this.$_app_aboutDialog = !1),
                (this.dialog = !0);
              var e = this;
              this.timer = setInterval(function () {
                e.timePassed += 1;
              }, 1e3);
            }
          },
          timePassed: function () {
            this.timeToReload <= this.timePassed && this.reload();
          },
        },
        beforeDestroy: function () {
          clearInterval(this.timer);
        },
        methods: {
          stopTimer: function () {
            clearInterval(this.timer), (this.dialog = !1);
          },
          reload: function () {
            this.stopTimer(), location.reload();
          },
        },
      },
      p = m,
      f = n("2877"),
      _ = n("6544"),
      v = n.n(_),
      g = n("8336"),
      x = n("b0af"),
      b = n("99d9"),
      y = n("12b2"),
      S = n("169a"),
      C = n("0e8f"),
      $ = n("132d"),
      T = n("a722"),
      w = n("3a2f"),
      D = Object(f["a"])(p, c, d, !1, null, null, null);
    D.options.__file = "TheConnectionStatus.vue";
    var k = D.exports;
    v()(D, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardActions: b["a"],
      VCardText: b["b"],
      VCardTitle: y["a"],
      VDialog: S["a"],
      VFlex: C["a"],
      VIcon: $["a"],
      VLayout: T["a"],
      VTooltip: w["a"],
    });
    var M = {
        components: { TheConnectionStatus: k },
        mixins: [u],
        computed: {
          disabled: function () {
            return (
              "running" === this.$_simulator_operationalStatus ||
              "paused" === this.$_simulator_operationalStatus
            );
          },
        },
      },
      E = M,
      P = n("9910"),
      R = n("71d9"),
      V = n("706c"),
      A = n("2a7f"),
      F = Object(f["a"])(E, r, l, !1, null, null, null);
    F.options.__file = "TheToolbar.vue";
    var O = F.exports;
    v()(F, {
      VSpacer: P["a"],
      VToolbar: R["a"],
      VToolbarSideIcon: V["a"],
      VToolbarTitle: A["a"],
    });
    var j = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "div",
          [
            n(
              "v-navigation-drawer",
              {
                attrs: {
                  fixed: "",
                  temporary: "",
                  "mobile-break-point": "false",
                  right: "",
                  app: "",
                },
                model: {
                  value: t.drawer,
                  callback: function (e) {
                    t.drawer = e;
                  },
                  expression: "drawer",
                },
              },
              [
                n(
                  "v-list",
                  { attrs: { dense: "" } },
                  t._l(t.items, function (e) {
                    return n(
                      "v-list-tile",
                      {
                        key: e.name,
                        attrs: { to: e.path, disabled: e.disabled },
                        on: {
                          click: function (t) {
                            return t.stopPropagation(), e.action(t);
                          },
                        },
                      },
                      [
                        n(
                          "v-list-tile-action",
                          [n("v-icon", [t._v(t._s(e.icon))])],
                          1
                        ),
                        n(
                          "v-list-tile-content",
                          [
                            n("v-list-tile-title", [
                              t._v(
                                "\n            " + t._s(e.name) + "\n          "
                              ),
                            ]),
                          ],
                          1
                        ),
                      ],
                      1
                    );
                  }),
                  1
                ),
              ],
              1
            ),
            n("TheShutdownDialog"),
          ],
          1
        );
      },
      I = [],
      L =
        (n("7f7f"),
        n("ac6a"),
        function () {
          var t = this,
            e = t.$createElement,
            n = t._self._c || e;
          return n(
            "v-dialog",
            {
              attrs: { persistent: "", "max-width": "400" },
              model: {
                value: t.$_app_shutdownDialog,
                callback: function (e) {
                  t.$_app_shutdownDialog = e;
                },
                expression: "$_app_shutdownDialog",
              },
            },
            [
              n(
                "v-card",
                [
                  n(
                    "v-layout",
                    [
                      n(
                        "v-card-title",
                        {
                          staticClass: "title",
                          attrs: { "primary-title": "" },
                        },
                        [t._v("\n        Shutdown the backend server?\n      ")]
                      ),
                    ],
                    1
                  ),
                  n(
                    "v-layout",
                    {
                      attrs: {
                        "align-center": "",
                        "justify-center": "",
                        row: "",
                      },
                    },
                    [
                      n(
                        "v-card-actions",
                        [
                          n(
                            "v-btn",
                            {
                              on: {
                                click: function (e) {
                                  e.stopPropagation(),
                                    (t.$_app_shutdownDialog = !1);
                                },
                              },
                            },
                            [t._v("\n          Cancel\n        ")]
                          ),
                          n(
                            "v-btn",
                            {
                              attrs: { color: "error" },
                              on: {
                                click: function (e) {
                                  e.stopPropagation(), t.shutdown();
                                },
                              },
                            },
                            [t._v("\n          Yes, do it now\n        ")]
                          ),
                        ],
                        1
                      ),
                    ],
                    1
                  ),
                ],
                1
              ),
            ],
            1
          );
        }),
      N = [],
      B = {
        mixins: [u],
        methods: {
          shutdown: function () {
            (this.$_app_shutdownDialog = !1),
              (this.$_app_autoReload = !1),
              this.$_eel.shutdown_backend()();
          },
        },
      },
      H = B,
      z = Object(f["a"])(H, L, N, !1, null, null, null);
    z.options.__file = "TheShutdownDialog.vue";
    var W = z.exports;
    v()(z, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardActions: b["a"],
      VCardTitle: y["a"],
      VDialog: S["a"],
      VLayout: T["a"],
    });
    var J = {
        components: { TheShutdownDialog: W },
        mixins: [u, h],
        data: function () {
          var t = this;
          return {
            items: [
              {
                name: "Dashboard",
                path: "/",
                icon: "home",
                action: function () {
                  t.close();
                },
                disabled: !0,
              },
              {
                name: "Settings",
                path: void 0,
                icon: "settings",
                action: function () {
                  (t.$_app_settingsDialog = !0), t.close();
                },
                disabled: !0,
              },
              {
                name: "Results",
                path: "/results",
                icon: "storage",
                action: function () {
                  t.close();
                },
                disabled: !0,
              },
              {
                name: "Trace Files",
                path: "/traces",
                icon: "folder",
                action: function () {
                  t.close();
                },
                disabled: !0,
              },
              {
                name: "Reload",
                path: void 0,
                icon: "refresh",
                action: function () {
                  location.reload();
                },
                disabled: !0,
              },
              {
                name: "Shutdown",
                path: void 0,
                icon: "power_settings_new",
                action: function () {
                  (t.$_app_shutdownDialog = !0), t.close();
                },
                disabled: !0,
              },
              {
                name: "About",
                path: void 0,
                icon: "info",
                action: function () {
                  (t.$_app_aboutDialog = !0), t.close();
                },
                disabled: !0,
              },
            ],
          };
        },
        computed: {
          drawer: {
            get: function () {
              return this.$store.getters["app/navigationDrawer"];
            },
            set: function (t) {
              t !== this.drawer &&
                this.$store.dispatch("app/toggleNavigationDrawer");
            },
          },
        },
        watch: {
          $_app_status: function (t) {
            this.items.forEach(function (e) {
              "Dashboard" === e.name || "About" === e.name
                ? (e.disabled = !1)
                : (e.disabled = "ready" !== t);
            });
          },
          $_simulator_operationalStatus: function (t) {
            this.items.forEach(function (e) {
              "Dashboard" === e.name || "About" === e.name
                ? (e.disabled = !1)
                : (e.disabled = "ready" !== t);
            });
          },
          $_simulator_connectionStatus: function (t) {
            this.items.forEach(function (e) {
              "disconnected" === t &&
                ("Reload" === e.name ? (e.disabled = !1) : (e.disabled = !0));
            });
          },
        },
        methods: {
          close: function () {
            this.drawer = !1;
          },
        },
      },
      U = J,
      K = n("8860"),
      X = n("ba95"),
      G = n("40fe"),
      q = n("5d23"),
      Y = n("f774"),
      Z = Object(f["a"])(U, j, I, !1, null, null, null);
    Z.options.__file = "TheNavigationDrawer.vue";
    var Q = Z.exports;
    v()(Z, {
      VIcon: $["a"],
      VList: K["a"],
      VListTile: X["a"],
      VListTileAction: G["a"],
      VListTileContent: q["a"],
      VListTileTitle: q["b"],
      VNavigationDrawer: Y["a"],
    });
    var tt = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-footer",
          { attrs: { height: "auto", app: "" } },
          [
            n(
              "v-flex",
              { attrs: { xs12: "" } },
              [
                n(
                  "v-container",
                  { attrs: { "pa-0": "", "ma-0": "", fluid: "" } },
                  [
                    n(
                      "v-layout",
                      { attrs: { "justify-center": "" } },
                      [
                        n(
                          "v-flex",
                          { attrs: { xs12: "" } },
                          [
                            n("v-progress-linear", {
                              staticClass: "ma-0",
                              attrs: {
                                indeterminate: t.running && t.progress < 5,
                                value: t.progress,
                                color: "pink lighten-3",
                                "background-color": "grey lighten-1",
                              },
                            }),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                    n(
                      "v-layout",
                      {
                        attrs: {
                          "align-center": "",
                          "justify-start": "",
                          row: "",
                          wrap: "",
                        },
                      },
                      [
                        n(
                          "v-btn",
                          {
                            attrs: {
                              icon: t.$vuetify.breakpoint.xsOnly,
                              disabled: t.disabled_play,
                            },
                            on: { click: t.onClickPlay },
                          },
                          [
                            t.$vuetify.breakpoint.smAndUp
                              ? n("span", [
                                  t._v(
                                    "\n            " +
                                      t._s(t._f("playString")(t.running)) +
                                      "\n          "
                                  ),
                                ])
                              : t._e(),
                            n("v-icon", [
                              t._v(t._s(t._f("playIcon")(t.running))),
                            ]),
                          ],
                          1
                        ),
                        n(
                          "v-btn",
                          {
                            attrs: {
                              icon: t.$vuetify.breakpoint.xsOnly,
                              disabled: t.disabled_stop,
                            },
                            on: { click: t.onClickStop },
                          },
                          [
                            t.$vuetify.breakpoint.smAndUp
                              ? n("span", [t._v("Stop")])
                              : t._e(),
                            n("v-icon", [t._v("stop")]),
                          ],
                          1
                        ),
                        n("v-flex", { attrs: { xs2: "", sm3: "" } }, [
                          n("span", [
                            n("span", [
                              t._v(
                                t._s(t._f("minutesString")(t.currentMinutes))
                              ),
                            ]),
                            t.$vuetify.breakpoint.smAndUp
                              ? n("span", [
                                  t._v(
                                    "\n              / " +
                                      t._s(
                                        t._f("minutesString")(t.endMinutes)
                                      ) +
                                      "\n            "
                                  ),
                                ])
                              : t._e(),
                          ]),
                        ]),
                        n("v-spacer"),
                        n(
                          "v-tooltip",
                          { attrs: { top: "" } },
                          [
                            n(
                              "v-chip",
                              {
                                attrs: { slot: "activator" },
                                slot: "activator",
                              },
                              [t._v("\n            " + t._s(t.sfClass))]
                            ),
                            n("span", [t._v("Scheduling Function")]),
                          ],
                          1
                        ),
                        n(
                          "v-tooltip",
                          { attrs: { top: "" } },
                          [
                            t.$vuetify.breakpoint.smAndUp
                              ? n(
                                  "v-chip",
                                  {
                                    attrs: { slot: "activator" },
                                    slot: "activator",
                                  },
                                  [
                                    t._v(
                                      "\n            " +
                                        t._s(t.connClass) +
                                        "\n          "
                                    ),
                                  ]
                                )
                              : t._e(),
                            n("span", [t._v("Connectivity Class")]),
                          ],
                          1
                        ),
                        n(
                          "v-tooltip",
                          { attrs: { top: "" } },
                          [
                            t.$vuetify.breakpoint.smAndUp &&
                            "K7" === t.connClass
                              ? n(
                                  "v-chip",
                                  {
                                    attrs: { slot: "activator" },
                                    slot: "activator",
                                  },
                                  [
                                    t._v(
                                      "\n            " +
                                        t._s(t.connTrace) +
                                        "\n          "
                                    ),
                                  ]
                                )
                              : t._e(),
                            n("span", [t._v("Connectivity Class")]),
                          ],
                          1
                        ),
                        n(
                          "v-tooltip",
                          { attrs: { top: "" } },
                          [
                            t.$vuetify.breakpoint.smAndUp
                              ? n(
                                  "v-chip",
                                  {
                                    attrs: { slot: "activator" },
                                    slot: "activator",
                                  },
                                  [
                                    t._v(
                                      "\n            " +
                                        t._s(t.numMotes) +
                                        "\n          "
                                    ),
                                  ]
                                )
                              : t._e(),
                            n("span", [t._v("Number of Motes")]),
                          ],
                          1
                        ),
                        n(
                          "v-tooltip",
                          { attrs: { top: "" } },
                          [
                            n(
                              "v-btn",
                              {
                                attrs: {
                                  slot: "activator",
                                  icon: "",
                                  disabled: !t.settingsButtonEnabled,
                                },
                                on: {
                                  click: function (e) {
                                    e.stopPropagation(),
                                      (t.$_app_settingsDialog = !0);
                                  },
                                },
                                slot: "activator",
                              },
                              [n("v-icon", [t._v("settings")])],
                              1
                            ),
                            n("span", [
                              t._v(
                                "\n            Change the simulator settings; you cannot do that while a\n            simulation is underway\n          "
                              ),
                            ]),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      et = [],
      nt =
        (n("6b54"),
        n("28a5"),
        n("4917"),
        n("a481"),
        {
          computed: {
            $_simulation_elapsedMinutes: function () {
              return this.$store.getters["simulation/elapsedMinutes"];
            },
            $_simulation_lastlastAppPacketEvent: function () {
              return this.$store.getters["simulation/lastAppPacketEvent"];
            },
            $_simulation_lastTschCellAllocationEvent: function () {
              return this.$store.getters[
                "simulation/lastTschCellAllocationEvent"
              ];
            },
            $_simulation_lastPacketDropEvent: function () {
              return this.$store.getters["simulation/lastPacketDropEvent"];
            },
            $_simulation_lastRplParentChangeEvent: function () {
              return this.$store.getters["simulation/lastRplParentChangeEvent"];
            },
            $_simulation_lastTschSyncEvent: function () {
              return this.$store.getters["simulation/lastTschSyncEvent"];
            },
            $_simulation_lastSecJoinEvent: function () {
              return this.$store.getters["simulation/lastSecJoinEvent"];
            },
          },
        }),
      at = {
        mixins: [u, h, nt],
        computed: {
          running: function () {
            return "running" === this.$_simulator_operationalStatus;
          },
          disabled_play: function () {
            return "ready" !== this.$_app_status || "/" !== this.$route.path;
          },
          disabled_stop: function () {
            return (
              "ready" !== this.$_app_status ||
              "ready" === this.$_simulator_operationalStatus ||
              "/" !== this.$route.path
            );
          },
          ready: function () {
            return "ready" === this.$_simulator_operationalStatus;
          },
          sfClass: function () {
            return null === this.$_simulator_runningSettings
              ? "unknown"
              : this.$_simulator_runningSettings.sf_class;
          },
          connClass: function () {
            return null === this.$_simulator_runningSettings
              ? "unknown"
              : this.$_simulator_runningSettings.conn_class;
          },
          connTrace: function () {
            if (
              null === this.$_simulator_runningSettings ||
              void 0 === this.$_simulator_runningSettings.conn_trace ||
              null === this.$_simulator_runningSettings.conn_trace
            )
              return "unknown";
            var t = this.$_simulator_runningSettings.conn_trace,
              e = t.replace(/^.*[\/]/, "");
            return e.match(/.k7.gz$/) ? e.split(".").slice(0, -2).join(".") : e;
          },
          numMotes: function () {
            return null === this.$_simulator_runningSettings
              ? "unknown"
              : this.$_simulator_runningSettings.exec_numMotes;
          },
          currentMinutes: function () {
            var t;
            return (
              (t =
                null === this.$_simulation_elapsedMinutes
                  ? 0
                  : this.$_simulation_elapsedMinutes),
              t
            );
          },
          endMinutes: function () {
            var t,
              e = this.$_simulator_runningSettings;
            return (
              (t =
                null === e
                  ? 0
                  : Math.floor(
                      (e.exec_numSlotframesPerRun *
                        e.tsch_slotframeLength *
                        e.tsch_slotDuration) /
                        60
                    )),
              t
            );
          },
          progress: function () {
            return 0 === this.endMinutes
              ? 0
              : (this.currentMinutes / this.endMinutes) * 100;
          },
          settingsButtonEnabled: function () {
            return (
              "ready" === this.$_simulator_operationalStatus &&
              null !== this.$_simulator_runningSettings &&
              this.$_simulator_availableSFs.length > 0 &&
              this.$_simulator_availableConnectivities.length > 0 &&
              "ready" === this.$_app_status
            );
          },
        },
        filters: {
          playString: function (t) {
            return !0 === t ? "Pause" : "Play";
          },
          playIcon: function (t) {
            return !0 === t ? "pause" : "play_arrow";
          },
          minutesString: function (t) {
            var e = ("0" + Math.floor(t / 60).toString()).slice(-2),
              n = ("0" + (t % 60).toString()).slice(-2);
            return e + "h" + n + "m";
          },
        },
        methods: {
          onClickPlay: function () {
            "running" === this.$_simulator_operationalStatus
              ? this.$_simulator_pause()
              : "paused" === this.$_simulator_operationalStatus
              ? this.$_simulator_resume()
              : this.$_simulator_start();
          },
          onClickStop: function () {
            this.$_simulator_abort();
          },
        },
      },
      st = at,
      it = n("cc20"),
      ot = n("a523"),
      rt = n("553a"),
      lt = n("8e36"),
      ut = Object(f["a"])(st, tt, et, !1, null, null, null);
    ut.options.__file = "TheFooter.vue";
    var ct = ut.exports;
    v()(ut, {
      VBtn: g["a"],
      VChip: it["a"],
      VContainer: ot["a"],
      VFlex: C["a"],
      VFooter: rt["a"],
      VIcon: $["a"],
      VLayout: T["a"],
      VProgressLinear: lt["a"],
      VSpacer: P["a"],
      VTooltip: w["a"],
    });
    var dt = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-dialog",
          {
            attrs: { persistent: "", width: "400px" },
            model: {
              value: t.$_app_aboutDialog,
              callback: function (e) {
                t.$_app_aboutDialog = e;
              },
              expression: "$_app_aboutDialog",
            },
          },
          [
            n(
              "v-card",
              [
                n(
                  "v-flex",
                  { attrs: { xs12: "" } },
                  [
                    n(
                      "v-container",
                      { attrs: { "grid-list-md": "" } },
                      [
                        n(
                          "v-layout",
                          { attrs: { row: "", wrap: "", "fill-height": "" } },
                          t._l(t.items, function (e) {
                            return n(
                              "v-flex",
                              { key: e.name, attrs: { xs12: "" } },
                              [
                                n(
                                  "v-card",
                                  { attrs: { flat: "" } },
                                  [
                                    n("v-card-title", [t._v(t._s(e.title))]),
                                    n("v-divider"),
                                    n(
                                      "v-card-text",
                                      [
                                        n(
                                          "v-layout",
                                          [
                                            n(
                                              "v-flex",
                                              { attrs: { xs2: "" } },
                                              [t._v("Path")]
                                            ),
                                            n(
                                              "v-flex",
                                              {
                                                attrs: {
                                                  xs10: "",
                                                  "text-xs-right": "",
                                                },
                                              },
                                              [
                                                t._v(
                                                  "\n                    " +
                                                    t._s(e.repo) +
                                                    "\n                  "
                                                ),
                                              ]
                                            ),
                                          ],
                                          1
                                        ),
                                        n(
                                          "v-layout",
                                          [
                                            n(
                                              "v-flex",
                                              { attrs: { xs2: "" } },
                                              [t._v("Branch")]
                                            ),
                                            n(
                                              "v-flex",
                                              {
                                                attrs: {
                                                  xs10: "",
                                                  "text-xs-right": "",
                                                },
                                              },
                                              [
                                                t._v(
                                                  "\n                    " +
                                                    t._s(e.branch) +
                                                    "\n                  "
                                                ),
                                              ]
                                            ),
                                          ],
                                          1
                                        ),
                                        n(
                                          "v-layout",
                                          [
                                            n(
                                              "v-flex",
                                              { attrs: { xs2: "" } },
                                              [t._v("Commit")]
                                            ),
                                            n(
                                              "v-flex",
                                              {
                                                attrs: {
                                                  xs10: "",
                                                  "text-xs-right": "",
                                                },
                                              },
                                              [
                                                t._v(
                                                  "\n                    " +
                                                    t._s(e.short_hash) +
                                                    "\n                  "
                                                ),
                                              ]
                                            ),
                                          ],
                                          1
                                        ),
                                      ],
                                      1
                                    ),
                                  ],
                                  1
                                ),
                              ],
                              1
                            );
                          }),
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
                n(
                  "v-card-actions",
                  [
                    n("v-spacer"),
                    n(
                      "v-btn",
                      {
                        on: {
                          click: function (e) {
                            return e.stopPropagation(), t.close(e);
                          },
                        },
                      },
                      [t._v("\n        OK\n      ")]
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      ht = [],
      mt =
        (n("f751"),
        n("456d"),
        {
          mixins: [u, h],
          computed: {
            items: function () {
              var t,
                e = this.$_simulator_gitInfo;
              return (
                (t =
                  null === e
                    ? [
                        {
                          name: "webapp",
                          title: "6TiSCH Simulator WebApp",
                          repo: "unknown",
                          branch: "unknown",
                          commit: "unknown",
                        },
                        {
                          name: "simulator",
                          title: "6TiSCH Simulator",
                          repo: "unknown",
                          branch: "unknown",
                          commit: "unknown",
                        },
                      ]
                    : Object.keys(e).map(function (t) {
                        var n;
                        return (
                          (n =
                            "webapp" === t
                              ? "6TiSCH Simulator WebApp"
                              : "simulator" === t
                              ? "6TiSCH Simulator"
                              : t),
                          Object.assign({ name: t, title: n }, e[t])
                        );
                      })),
                t
              );
            },
          },
          methods: {
            close: function () {
              this.$_app_aboutDialog = !1;
            },
          },
        }),
      pt = mt,
      ft = n("ce7e"),
      _t = Object(f["a"])(pt, dt, ht, !1, null, null, null);
    _t.options.__file = "TheAboutDialog.vue";
    var vt = _t.exports;
    v()(_t, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardActions: b["a"],
      VCardText: b["b"],
      VCardTitle: y["a"],
      VContainer: ot["a"],
      VDialog: S["a"],
      VDivider: ft["a"],
      VFlex: C["a"],
      VLayout: T["a"],
      VSpacer: P["a"],
    });
    var gt = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "div",
          [
            n(
              "v-dialog",
              {
                attrs: { persistent: "", width: "500px" },
                model: {
                  value: t.$_app_settingsDialog,
                  callback: function (e) {
                    t.$_app_settingsDialog = e;
                  },
                  expression: "$_app_settingsDialog",
                },
              },
              [
                n(
                  "v-tabs",
                  {
                    attrs: { "fixed-tabs": "", "icons-and-text": "" },
                    model: {
                      value: t.active,
                      callback: function (e) {
                        t.active = e;
                      },
                      expression: "active",
                    },
                  },
                  [
                    t._l(t.items, function (e) {
                      return n(
                        "v-tab",
                        { key: e.name },
                        [
                          t._v("\n        " + t._s(e.name) + "\n        "),
                          n("v-icon", [t._v(t._s(e.icon))]),
                        ],
                        1
                      );
                    }),
                    n(
                      "v-tab-item",
                      { key: t.items[0].name },
                      [n("TheSimpleConfiguratorTab")],
                      1
                    ),
                    n(
                      "v-tab-item",
                      { key: t.items[1].name },
                      [n("TheConfigFileUploadTab")],
                      1
                    ),
                    n(
                      "v-tab-item",
                      { key: t.items[2].name },
                      [n("TheConfigFileDownloadTab")],
                      1
                    ),
                  ],
                  2
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      xt = [],
      bt = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-card",
          [
            n(
              "v-card-text",
              [
                n(
                  "v-list",
                  t._l(t.items, function (e) {
                    return n(
                      "v-list-tile",
                      { key: e.name, attrs: { disabled: e.disabled } },
                      [
                        n("v-list-tile-title", [t._v(t._s(e.title))]),
                        "number" === e.type
                          ? n("v-text-field", {
                              attrs: {
                                "single-line": "",
                                type: "number",
                                disabled: e.disabled,
                                rules: e.rules,
                                reverse: "",
                              },
                              model: {
                                value: e.model,
                                callback: function (n) {
                                  t.$set(e, "model", t._n(n));
                                },
                                expression: "item.model",
                              },
                            })
                          : n("v-select", {
                              staticClass: "text-xs-right",
                              attrs: {
                                items: e.selectItems,
                                disabled: e.disabled,
                                reverse: "",
                              },
                              model: {
                                value: e.model,
                                callback: function (n) {
                                  t.$set(e, "model", n);
                                },
                                expression: "item.model",
                              },
                            }),
                      ],
                      1
                    );
                  }),
                  1
                ),
              ],
              1
            ),
            n(
              "v-card-actions",
              [
                n("v-spacer"),
                n(
                  "v-btn",
                  {
                    attrs: { color: "error" },
                    on: {
                      click: function (e) {
                        return e.stopPropagation(), t.resetSettings(e);
                      },
                    },
                  },
                  [t._v("Reset")]
                ),
                n(
                  "v-btn",
                  {
                    attrs: { color: "" },
                    on: {
                      click: function (e) {
                        return e.stopPropagation(), t.cancelSettings(e);
                      },
                    },
                  },
                  [t._v("Cancel")]
                ),
                n(
                  "v-btn",
                  {
                    attrs: { color: "primary", disabled: t.disabledSave },
                    on: {
                      click: function (e) {
                        return e.stopPropagation(), t.saveSettings(e);
                      },
                    },
                  },
                  [t._v("\n      Save\n    ")]
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      yt = [],
      St =
        (n("7514"),
        {
          mixins: [u, h],
          data: function () {
            var t = this;
            return {
              items: [
                {
                  name: "sf_class",
                  title: "Scheduling Function",
                  model: "",
                  type: "select",
                  selectItems: [],
                  disabled: !1,
                  ready: !1,
                },
                {
                  name: "conn_class",
                  title: "Connectivity Class",
                  model: "",
                  type: "select",
                  selectItems: [],
                  disabled: !1,
                  ready: !1,
                },
                {
                  name: "conn_trace",
                  title: "Trace File",
                  model: "",
                  type: "select",
                  selectItems: [],
                  disabled: !0,
                  ready: !1,
                },
                {
                  name: "exec_numMotes",
                  title: "Number of Motes",
                  model: "",
                  modelMax: 1 / 0,
                  type: "number",
                  disabled: !1,
                  rules: [
                    function (t) {
                      return "" !== t || "input is required";
                    },
                    function (t) {
                      return t > 0 || "input must be larger than 0";
                    },
                  ],
                  ready: !1,
                },
                {
                  name: "exec_randomSeed",
                  title: "Random Seed",
                  model: "",
                  type: "number",
                  disabled: !1,
                  rules: [
                    function (t) {
                      return "" !== t || "input is required";
                    },
                  ],
                  ready: !1,
                },
                {
                  name: "exec_numMinutes",
                  title: "Simulation Time (minutes)",
                  model: "",
                  modelMax: 1 / 0,
                  type: "number",
                  disabled: !1,
                  rules: [
                    function (t) {
                      return "" !== t || "input is required";
                    },
                    function (t) {
                      return t > 0 || "input must be larger than 0";
                    },
                    function (e) {
                      return (
                        e <= t.exec_numMinutes.modelMax ||
                        "input must not exceed " + t.exec_numMinutes.modelMax
                      );
                    },
                  ],
                  ready: !1,
                },
              ],
              exec_numMotes_backup: 0,
            };
          },
          computed: {
            disabledSave: function () {
              return !this.items.every(function (t) {
                return !0 === t.ready;
              });
            },
            sf_class: function () {
              return this.items[0];
            },
            conn_class: function () {
              return this.items[1];
            },
            conn_trace: function () {
              return this.items[2];
            },
            exec_numMotes: function () {
              return this.items[3];
            },
            exec_randomSeed: function () {
              return this.items[4];
            },
            exec_numMinutes: function () {
              return this.items[5];
            },
            exec_numSlotframesPerRun: {
              get: function () {
                if (null === this.$_simulator_runningSettings) return 0;
                var t = this.$_simulator_runningSettings;
                return Math.ceil(
                  (60 * this.exec_numMinutes_model) /
                    t.tsch_slotDuration /
                    t.tsch_slotframeLength
                );
              },
              set: function (t) {
                var e = this.$_simulator_runningSettings,
                  n = Math.floor(
                    (t * e.tsch_slotframeLength * e.tsch_slotDuration) / 60
                  );
                return (this.exec_numMinutes_model = n);
              },
            },
            sf_class_model: {
              get: function () {
                return this.sf_class.model;
              },
              set: function (t) {
                this.sf_class.model = t;
              },
            },
            conn_class_model: {
              get: function () {
                return this.conn_class.model;
              },
              set: function (t) {
                this.conn_class.model = t;
              },
            },
            conn_trace_model: {
              get: function () {
                return this.conn_trace.model;
              },
              set: function (t) {
                this.conn_trace.model = t;
              },
            },
            exec_numMotes_model: {
              get: function () {
                return this.exec_numMotes.model;
              },
              set: function (t) {
                this.exec_numMotes.model = t;
              },
            },
            exec_randomSeed_model: {
              get: function () {
                return this.exec_randomSeed.model;
              },
              set: function (t) {
                this.exec_randomSeed.model = t;
              },
            },
            exec_numMinutes_model: {
              get: function () {
                return this.exec_numMinutes.model;
              },
              set: function (t) {
                this.exec_numMinutes.model = t;
              },
            },
            newSettings: function () {
              return Object.assign(this.$_simulator_runningSettings, {
                sf_class: this.sf_class_model,
                conn_class: this.conn_class_model,
                conn_trace: this.conn_trace_model,
                exec_numMotes: this.exec_numMotes_model,
                exec_randomSeed: this.exec_randomSeed_model,
                exec_numSlotframesPerRun: this.exec_numSlotframesPerRun,
              });
            },
          },
          watch: {
            sf_class_model: function (t) {
              this.sf_class.ready = "" !== t;
            },
            conn_class_model: function (t) {
              if (((this.conn_class.ready = "" !== t), "K7" === t)) {
                this.conn_trace.disabled = !1;
                var e = this.$_simulator_runningSettings.conn_trace;
                (this.conn_trace_model =
                  null === e ? this.conn_trace.selectItems[0].value : e),
                  (this.exec_numMotes.disabled = !0),
                  "" === this.conn_trace_model
                    ? (this.conn_trace.ready = !1)
                    : (this.conn_trace.ready = !0);
              } else
                (this.conn_trace.disabled = !0),
                  (this.conn_trace_model = ""),
                  (this.exec_numMotes.disabled = !1),
                  (this.conn_trace.ready = !0);
            },
            conn_trace_model: function (t) {
              "K7" !== this.conn_class_model || "" !== t
                ? (this.conn_trace.ready = !0)
                : (this.conn_trace.ready = !1),
                "" === t
                  ? this.exec_numMotes_backup > 0 &&
                    ((this.exec_numMotes_model = this.exec_numMotes_backup),
                    (this.exec_numMotes_backup = 0),
                    (this.exec_numMinutes.modelMax = 1 / 0))
                  : ((this.exec_numMotes_backup = this.exec_numMotes_model),
                    (this.exec_numMotes_model = this.getNumMotesFromTrace(t)),
                    (this.exec_numMinutes.modelMax =
                      this.getMaxExecNumMinuetsFromTrace(t)));
            },
            exec_numMotes_model: function (t) {
              this.exec_numMotes.ready =
                t > 0 && t < this.exec_numMotes.modelMax;
            },
            exec_randomSeed_model: function (t) {
              this.exec_randomSeed.ready = "" !== t;
            },
            exec_numMinutes_model: function (t) {
              this.exec_numMinutes.ready =
                t > 0 && t < this.exec_numMinutes.modelMax;
            },
            $_simulator_runningSettings: function (t) {
              null !== t &&
                ((this.sf_class_model = t.sf_class),
                (this.conn_class_model = t.conn_class),
                (this.exec_numMotes_model = t.exec_numMotes),
                (this.exec_randomSeed_model = t.exec_randomSeed),
                (this.exec_numSlotframesPerRun = t.exec_numSlotframesPerRun));
            },
            $_simulator_availableSFs: function (t) {
              this.sf_class.selectItems = null === t ? [] : t;
            },
            $_simulator_availableConnectivities: function (t) {
              this.conn_class.selectItems = null === t ? [] : t;
            },
            $_simulator_availableTraceFiles: function (t) {
              if (null === t) this.conn_trace.selectItems = [];
              else if (
                ((this.conn_trace.selectItems = t.map(function (t) {
                  var e;
                  return (
                    (e = t.file_name.match(/.k7.gz$/)
                      ? t.file_name.split(".").slice(0, -2).join(".")
                      : t.file_name),
                    { text: e, value: t.file_path }
                  );
                })),
                "" !== this.conn_trace_model)
              ) {
                var e = this.conn_trace_model;
                (this.exec_numMotes_model = this.getNumMotesFromTrace(e)),
                  (this.exec_numMinutes.modelMax =
                    this.getMaxExecNumMinuetsFromTrace(e));
              }
            },
          },
          methods: {
            saveSettings: function () {
              (this.$_simulator_runningSettings = this.newSettings),
                (this.$_app_settingsDialog = !1);
            },
            resetSettings: function () {
              if (null === this.$_simulator_runningSettings)
                (this.sf_class_model = ""),
                  (this.conn_class_model = ""),
                  (this.exec_numMotes_model = ""),
                  (this.exec_randomSeed_model = ""),
                  (this.exec_numSlotframesPerRun = "");
              else {
                var t = this.$_simulator_runningSettings;
                (this.sf_class_model = t.sf_class),
                  (this.conn_class_model = t.conn_class),
                  (this.exec_numMotes_model = t.exec_numMotes),
                  (this.exec_randomSeed_model = t.exec_randomSeed),
                  (this.exec_numSlotframesPerRun = t.exec_numSlotframesPerRun);
              }
            },
            cancelSettings: function () {
              this.resetSettings(), (this.$_app_settingsDialog = !1);
            },
            getNumMotesFromTrace: function (t) {
              var e = this.getTrace(t);
              return void 0 === e || null === e.config
                ? 0
                : e.config.node_count;
            },
            getMaxExecNumMinuetsFromTrace: function (t) {
              var e = this.getTrace(t);
              return void 0 === e ? 1 / 0 : e.config.maxDuration;
            },
            getTrace: function (t) {
              return this.$_simulator_availableTraceFiles.find(function (e) {
                return e.file_path === t;
              });
            },
            checkMaxExecNumMinutes: function () {},
          },
        }),
      Ct = St,
      $t = n("b56d"),
      Tt = n("2677"),
      wt = Object(f["a"])(Ct, bt, yt, !1, null, null, null);
    wt.options.__file = "TheSimpleConfiguratorTab.vue";
    var Dt = wt.exports;
    v()(wt, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardActions: b["a"],
      VCardText: b["b"],
      VList: K["a"],
      VListTile: X["a"],
      VListTileTitle: q["b"],
      VSelect: $t["a"],
      VSpacer: P["a"],
      VTextField: Tt["a"],
    });
    var kt = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "div",
          [
            n(
              "v-card",
              [
                n(
                  "v-layout",
                  { attrs: { "justify-center": "", row: "" } },
                  [
                    n(
                      "v-flex",
                      { attrs: { xs10: "" } },
                      [
                        n(
                          "v-alert",
                          {
                            attrs: {
                              value: t.alert,
                              outline: "",
                              type: "error",
                            },
                          },
                          [
                            t._v(
                              "\n          Upload Failed: " +
                                t._s(t.errorMessage) +
                                "\n        "
                            ),
                          ]
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
                n("v-card-text", { staticClass: "title text-xs-center" }, [
                  t._v(" " + t._s(t.selectedFileName) + "\n    "),
                ]),
                n(
                  "v-layout",
                  {
                    attrs: {
                      "align-center": "",
                      "justify-center": "",
                      row: "",
                      wrap: "",
                      "fill-height": "",
                    },
                  },
                  [
                    n("v-btn", { on: { click: t.selectFile } }, [
                      t._v(
                        "\n        Select config.json on your computer\n      "
                      ),
                    ]),
                    n("input", {
                      ref: "fileSelector",
                      attrs: {
                        type: "file",
                        multiple: "false",
                        accept: "application/json,.json",
                        hidden: "",
                      },
                      on: {
                        change: function (e) {
                          t.changeSelectedFile();
                        },
                      },
                    }),
                    n(
                      "v-btn",
                      {
                        on: {
                          click: function (e) {
                            t.$_app_settingsDialog = !1;
                          },
                        },
                      },
                      [t._v("Cancel")]
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
            n(
              "v-dialog",
              {
                attrs: { persistent: "", width: "300" },
                model: {
                  value: t.progress,
                  callback: function (e) {
                    t.progress = e;
                  },
                  expression: "progress",
                },
              },
              [
                n(
                  "v-card",
                  { attrs: { color: "primary", dark: "" } },
                  [
                    n(
                      "v-card-text",
                      [
                        t._v("\n        Please stand by\n        "),
                        n("v-progress-linear", {
                          staticClass: "mb-0",
                          attrs: { indeterminate: "", color: "white" },
                        }),
                      ],
                      1
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      Mt = [],
      Et = {
        mixins: [u, h],
        data: function () {
          return { progress: !1, errorMessage: null, selectedFile: null };
        },
        computed: {
          selector: function () {
            return this.$refs.fileSelector.value;
          },
          alert: function () {
            return null !== this.errorMessage;
          },
          selectedFileName: function () {
            return null === this.selectedFile
              ? "No file is selected"
              : this.selectedFile.name;
          },
        },
        methods: {
          selectFile: function () {
            this.$refs.fileSelector.click();
          },
          changeSelectedFile: function () {
            1 === this.$refs.fileSelector.files.length &&
              ((this.errorMessage = null),
              (this.selectedFile = this.$refs.fileSelector.files[0]),
              this.uploadConfigFile());
          },
          uploadConfigFile: function () {
            var t = this,
              e = new FileReader();
            (this.progress = !0),
              (e.onload = function (e) {
                (t.progress = !1),
                  t.$_eel.put_default_config(e.target.result)(function (e) {
                    (t.progress = !1),
                      null === e.config
                        ? (t.errorMessage = e.message)
                        : ((t.$_simulator_defaultSettings = e.config.settings),
                          (t.$_app_settingsDialog = !1)),
                      (t.selectedFile = null);
                  });
              }),
              e.readAsText(this.selectedFile);
          },
        },
      },
      Pt = Et,
      Rt = n("0798"),
      Vt = Object(f["a"])(Pt, kt, Mt, !1, null, null, null);
    Vt.options.__file = "TheConfigFileUploadTab.vue";
    var At = Vt.exports;
    v()(Vt, {
      VAlert: Rt["a"],
      VBtn: g["a"],
      VCard: x["a"],
      VCardText: b["b"],
      VDialog: S["a"],
      VFlex: C["a"],
      VLayout: T["a"],
      VProgressLinear: lt["a"],
    });
    var Ft = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-card",
          [
            n("v-card-text", { staticClass: "title text-xs-center" }, [
              t._v("Click to Download config.json"),
            ]),
            n(
              "v-layout",
              { attrs: { "justify-center": "", row: "", "fill-height": "" } },
              [
                n(
                  "v-btn",
                  {
                    attrs: { icon: "", href: "/config.json" },
                    on: {
                      click: function (e) {
                        t.$_app_settingsDialog = !1;
                      },
                    },
                  },
                  [n("v-icon", { attrs: { large: "" } }, [t._v("save_alt")])],
                  1
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      Ot = [],
      jt = { mixins: [u] },
      It = jt,
      Lt = Object(f["a"])(It, Ft, Ot, !1, null, null, null);
    Lt.options.__file = "TheConfigFileDownloadTab.vue";
    var Nt = Lt.exports;
    v()(Lt, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardText: b["b"],
      VIcon: $["a"],
      VLayout: T["a"],
    });
    var Bt = {
        components: {
          TheSimpleConfiguratorTab: Dt,
          TheConfigFileUploadTab: At,
          TheConfigFileDownloadTab: Nt,
        },
        mixins: [u],
        data: function () {
          return {
            active: null,
            items: [
              { name: "Edit", icon: "edit" },
              { name: "Upload", icon: "cloud_upload" },
              { name: "Download", icon: "cloud_download" },
            ],
          };
        },
      },
      Ht = Bt,
      zt = n("71a3"),
      Wt = n("c671"),
      Jt = n("fe57"),
      Ut = Object(f["a"])(Ht, gt, xt, !1, null, null, null);
    Ut.options.__file = "TheSettingsDialog.vue";
    var Kt = Ut.exports;
    v()(Ut, {
      VDialog: S["a"],
      VIcon: $["a"],
      VTab: zt["a"],
      VTabItem: Wt["a"],
      VTabs: Jt["a"],
    });
    var Xt = {
        name: "App",
        components: {
          TheToolbar: O,
          TheNavigationDrawer: Q,
          TheFooter: ct,
          TheAboutDialog: vt,
          TheSettingsDialog: Kt,
        },
        data: function () {
          return {};
        },
        created: function () {
          document.addEventListener("DOMContentLoaded", this.initializeEel);
        },
        methods: {
          initializeEel: function () {
            if (void 0 === this.$_eel) this.eelCloseHandler();
            else {
              var t = this.$_eel._websocket,
                e = this;
              t.readyState === t.CONNECTING
                ? t.addEventListener("open", function () {
                    e.eelOpenHandler();
                  })
                : t.readyState === t.OPEN && this.eelOpenHandler(),
                t.readyState === t.CLOSING || t.readyState === t.CLOSED
                  ? this.eelCloseHandler()
                  : t.addEventListener("close", function () {
                      e.eelCloseHandler();
                    });
            }
          },
          eelOpenHandler: function () {
            var t = this;
            this.$store.dispatch("simulator/connect"),
              this.$_eel.get_default_config()(function (e) {
                t.$store.dispatch("simulator/saveDefaultSettings", e.settings);
              }),
              this.$_eel.get_available_scheduling_functions()(function (e) {
                t.$store.dispatch("simulator/setAvailableSFs", e);
              }),
              this.$_eel.get_available_connectivities()(function (e) {
                t.$store.dispatch("simulator/setAvailableConnectivities", e);
              }),
              this.$_eel.get_available_trace_files()(function (e) {
                t.$store.dispatch("simulator/setAvailableTraceFiles", e);
              }),
              this.$_eel.get_git_info()(function (e) {
                t.$store.dispatch("simulator/setGitInfo", e);
              });
          },
          eelCloseHandler: function () {
            this.$store.dispatch("simulator/disconnect");
          },
        },
      },
      Gt = Xt,
      qt = n("7496"),
      Yt = n("549c"),
      Zt = Object(f["a"])(Gt, i, o, !1, null, null, null);
    Zt.options.__file = "App.vue";
    var Qt = Zt.exports;
    v()(Zt, { VApp: qt["a"], VContent: Yt["a"], VLayout: T["a"] });
    var te = n("2f62"),
      ee = {
        namespaced: !0,
        state: {
          navigationDrawer: !1,
          shutdownDialog: !1,
          autoReload: !0,
          aboutDialog: !1,
          settingsDialog: !1,
          status: "starting",
        },
        getters: {
          navigationDrawer: function (t) {
            return t.navigationDrawer;
          },
          shutdownDialog: function (t) {
            return t.shutdownDialog;
          },
          autoReload: function (t) {
            return t.autoReload;
          },
          aboutDialog: function (t) {
            return t.aboutDialog;
          },
          settingsDialog: function (t) {
            return t.settingsDialog;
          },
          status: function (t) {
            return t.status;
          },
        },
        mutations: {
          setNavigationDrawer: function (t, e) {
            t.navigationDrawer = e;
          },
          setShutdownDialog: function (t, e) {
            t.shutdownDialog = e;
          },
          setAutoReload: function (t, e) {
            t.autoReload = e;
          },
          setAboutDialog: function (t, e) {
            t.aboutDialog = e;
          },
          setConfigDialog: function (t, e) {
            t.settingsDialog = e;
          },
          setStatus: function (t, e) {
            t.status = e;
          },
        },
        actions: {
          toggleNavigationDrawer: function (t) {
            var e = t.state,
              n = t.commit;
            n("setNavigationDrawer", !e.navigationDrawer);
          },
          setShutdownDialog: function (t, e) {
            var n = t.commit;
            n("setShutdownDialog", e);
          },
          setAutoReload: function (t, e) {
            var n = t.commit;
            n("setAutoReload", e);
          },
          setAboutDialog: function (t, e) {
            var n = t.commit;
            n("setAboutDialog", e);
          },
          setConfigDialog: function (t, e) {
            var n = t.commit;
            n("setConfigDialog", e);
          },
          setStatus: function (t, e) {
            var n = t.commit;
            n("setStatus", e);
          },
        },
      },
      ne = {
        namespaced: !0,
        state: {
          defaultLogFilter: [
            "app.rx",
            "packet_dropped",
            "tsch.add_cell",
            "tsch.delete_cell",
            "tsch.synced",
            "tsch.desynced",
            "secjoin.joined",
            "secjoin.unjoined",
            "rpl.churn",
          ],
          connectionStatus: "disconnected",
          operationalStatus: null,
          defaultSettings: null,
          runningSettings: null,
          availableSFs: [],
          availableConnectivities: [],
          availableTraceFiles: [],
          gitInfo: null,
          crashReport: null,
        },
        getters: {
          connectionStatus: function (t) {
            return t.connectionStatus;
          },
          operationalStatus: function (t) {
            return t.operationalStatus;
          },
          defaultSettings: function (t) {
            return t.defaultSettings;
          },
          runningSettings: function (t) {
            return t.runningSettings;
          },
          availableSFs: function (t) {
            return t.availableSFs;
          },
          availableConnectivities: function (t) {
            return t.availableConnectivities;
          },
          availableTraceFiles: function (t) {
            return t.availableTraceFiles;
          },
          gitInfo: function (t) {
            return t.gitInfo;
          },
          crashReport: function (t) {
            return t.crashReport;
          },
        },
        mutations: {
          changeConnectionStatus: function (t, e) {
            t.connectionStatus = e;
          },
          changeOperationalStatus: function (t, e) {
            t.operationalStatus = e;
          },
          saveDefaultSettings: function (t, e) {
            t.defaultSettings = e;
          },
          updateRunningSettings: function (t, e) {
            t.runningSettings = Object.assign({}, e);
          },
          setAvailableSFs: function (t, e) {
            t.availableSFs = e;
          },
          setAvailableConnectivities: function (t, e) {
            t.availableConnectivities = e;
          },
          setAvailableTraceFiles: function (t, e) {
            t.availableTraceFiles = e;
          },
          setGitInfo: function (t, e) {
            t.gitInfo = e;
          },
          setCrashReport: function (t, e) {
            t.crashReport = e;
          },
        },
        actions: {
          disconnect: function (t) {
            var e = t.commit;
            e("changeConnectionStatus", "disconnected");
          },
          connect: function (t) {
            var e = t.commit;
            e("changeConnectionStatus", "connected");
          },
          ready: function (t) {
            var e = t.commit;
            e("changeOperationalStatus", "ready");
          },
          start: function (t, e) {
            var n = t.state,
              a = t.commit,
              s = t.dispatch;
            s("simulation/reset", null, { root: !0 }),
              a("changeOperationalStatus", "running"),
              e.start(
                n.runningSettings,
                n.defaultLogFilter
              )(function (t) {
                "success" === t.status ||
                  "aborted" === t.status ||
                  ("failure" === t.status && a("setCrashReport", t)),
                  s("ready");
              });
          },
          pause: function (t, e) {
            var n = t.commit;
            e.pause()(function () {
              n("changeOperationalStatus", "paused");
            });
          },
          resume: function (t, e) {
            var n = t.commit;
            e.resume()(function () {
              n("changeOperationalStatus", "running");
            });
          },
          abort: function (t, e) {
            var n = t.commit,
              a = t.dispatch;
            n("changeOperationalStatus", "aborted"),
              e.abort()(function () {
                a("simulation/reset", null, { root: !0 }), a("ready");
              });
          },
          saveDefaultSettings: function (t, e) {
            var n,
              a = t.commit,
              s = t.dispatch;
            if ((a("saveDefaultSettings", e), null === e)) n = null;
            else
              for (var i in ((n = Object.assign({}, e.regular)),
              "random" === e.regular.exec_randomSeed &&
                (n.exec_randomSeed = Math.floor(1e4 * Math.random())),
              e.combination))
                n[i] = e.combination[i][0];
            s("updateRunningSettings", n), s("ready");
          },
          updateRunningSettings: function (t, e) {
            var n = t.commit;
            n("updateRunningSettings", e);
          },
          setAvailableSFs: function (t, e) {
            var n = t.commit;
            n("setAvailableSFs", e);
          },
          setAvailableConnectivities: function (t, e) {
            var n = t.commit;
            n("setAvailableConnectivities", e);
          },
          setAvailableTraceFiles: function (t, e) {
            var n = t.commit;
            n(
              "setAvailableTraceFiles",
              e.map(function (t) {
                var e = new Date(t.config.start_date),
                  n = new Date(t.config.stop_date);
                return (
                  (t.config.maxDuration = Math.floor((n - e) / 1e3 / 60)), t
                );
              })
            );
          },
          setGitInfo: function (t, e) {
            var n = t.commit;
            n("setGitInfo", e);
          },
          clearCrashReport: function (t) {
            var e = t.commit;
            e("setCrashReport", null);
          },
        },
      },
      ae =
        (n("6c7b"),
        n("20d6"),
        {
          namespaced: !0,
          state: {
            lastRplParentChangeEvent: null,
            lastTschCellAllocationEvent: null,
            lastAppPacketEvent: null,
            lastPacketDropEvent: null,
            lastTschSyncEvent: null,
            lastSecJoinEvent: null,
            elapsedMinutes: null,
            motes: [],
          },
          getters: {
            lastRplParentChangeEvent: function (t) {
              return t.lastRplParentChangeEvent;
            },
            lastTschCellAllocationEvent: function (t) {
              return t.lastTschCellAllocationEvent;
            },
            lastAppPacketEvent: function (t) {
              return t.lastAppPacketEvent;
            },
            lastPacketDropEvent: function (t) {
              return t.lastPacketDropEvent;
            },
            lastTschSyncEvent: function (t) {
              return t.lastTschSyncEvent;
            },
            lastSecJoinEvent: function (t) {
              return t.lastSecJoinEvent;
            },
            elapsedMinutes: function (t) {
              return t.elapsedMinutes;
            },
          },
          mutations: {
            updateLastRplParentChangeEvent: function (t, e) {
              if (null !== e) {
                var n,
                  a = e._mote_id,
                  s = t.motes[a].rplParentId;
                (n =
                  null == e.preferredParent
                    ? null
                    : t.motes.findIndex(function (t) {
                        return t.eui64Addr === e.preferredParent;
                      })),
                  (t.motes[a].rplParentId = n),
                  (t.lastRplParentChangeEvent = {
                    childId: a,
                    oldParentId: s,
                    newParentId: n,
                  });
              }
            },
            updateLastTschCellAllocationEvent: function (t, e) {
              var n;
              null !== e &&
                ("tsch.add_cell" === e._type
                  ? (n = "add")
                  : "tsch.delete_cell" === e._type && (n = "delete"),
                (t.lastTschCellAllocationEvent = {
                  type: n,
                  moteId: e._mote_id,
                  slotFrameHandle: e.slotFrameHandle,
                  slotOffset: e.slotOffset,
                  channelOffset: e.channelOffset,
                  neighbor: e.neighbor,
                  cellOptions: e.cellOptions,
                }));
            },
            updateLastAppPacketEvent: function (t, e) {
              var n;
              null === e
                ? (t.lastAppPacketEvent = null)
                : ((n = "app.rx" === e._type ? "rx" : "drop"),
                  (t.lastAppPacketEvent = {
                    type: n,
                    asn: e._asn,
                    packet: e.packet,
                  }));
            },
            updateLastPacketDropEvent: function (t, e) {
              t.lastPacketDropEvent =
                null === e
                  ? null
                  : {
                      asn: e._asn,
                      moteId: e._mote_id,
                      packetType: e.packet.type,
                      dropReason: e.reason,
                    };
            },
            updateLastTschSyncEvent: function (t, e) {
              t.lastTschSyncEvent = e;
            },
            updateLastSecJoinEvent: function (t, e) {
              t.lastSecJoinEvent = e;
            },
            updateElapsedMinutes: function (t, e) {
              t.elapsedMinutes = e;
            },
            updateMoteAddress: function (t, e) {
              null !== e &&
                ("mac.add_addr" === e._type && "eui64" === e.type
                  ? (t.motes[e._mote_id].eui64Addr = e.addr)
                  : "ipv6.add_addr" === e._type &&
                    ("global" === e.type
                      ? (t.motes[e._mote_id].ipv6GlobalAddr = e.addr)
                      : (t.motes[e._mote_id].ipv6LinkLocalAddr = e.addr)));
            },
            initializeMotes: function (t, e) {
              t.motes = Array(e)
                .fill()
                .map(function () {
                  return {
                    eui64Addr: null,
                    ipv6LinkLocalAddr: null,
                    ipv6GlobalAddr: null,
                    rplParentId: null,
                  };
                });
            },
          },
          actions: {
            put: function (t, e) {
              var n = t.commit;
              "_backend.tick.minute" === e._type
                ? n("updateElapsedMinutes", e.currentValue)
                : "app.rx" === e._type
                ? n("updateLastAppPacketEvent", e)
                : "packet_dropped" === e._type
                ? ("DATA" === e.packet.type && n("updateLastAppPacketEvent", e),
                  n("updateLastPacketDropEvent", e))
                : "mac.add_addr" === e._type && "eui64" === e.type
                ? n("updateMoteAddress", e)
                : "tsch.add_cell" === e._type || "tsch.delete_cell" === e._type
                ? n("updateLastTschCellAllocationEvent", e)
                : "rpl.churn" === e._type
                ? n("updateLastRplParentChangeEvent", e)
                : "tsch.synced" === e._type || "tsch.desynced" === e._type
                ? n("updateLastTschSyncEvent", e)
                : ("secjoin.joined" !== e._type &&
                    "secjoin.unjoined" !== e._type) ||
                  n("updateLastSecJoinEvent", e);
            },
            reset: function (t) {
              var e = t.commit,
                n = t.rootGetters;
              e("updateElapsedMinutes", 0),
                e("updateLastAppPacketEvent", null),
                e("updateLastPacketDropEvent", null),
                e("updateLastTschCellAllocationEvent", null),
                e("updateLastRplParentChangeEvent", null),
                e(
                  "initializeMotes",
                  n["simulator/runningSettings"].exec_numMotes
                );
            },
          },
        });
    a["a"].use(te["a"]);
    var se = new te["a"].Store({
        modules: { app: ee, simulator: ne, simulation: ae },
      }),
      ie = n("8c4f"),
      oe = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs12: "" } },
          [
            n(
              "v-container",
              { attrs: { "grid-list-sm": "" } },
              [
                n(
                  "v-layout",
                  { attrs: { row: "", wrap: "" } },
                  [
                    n("TheOverallE2ePdrCard"),
                    n("TheOverallE2eLatencyCard"),
                    n("TheSlotframeOccupancyCard"),
                    n("TheMoteStatsCard"),
                  ],
                  1
                ),
                n(
                  "v-layout",
                  { attrs: { row: "", wrap: "" } },
                  [
                    n("TheE2ePdrChartCard"),
                    n("TheE2eLatencyChartCard"),
                    n("ThePacketDropChartCard"),
                  ],
                  1
                ),
                n(
                  "v-layout",
                  { attrs: { row: "", wrap: "" } },
                  [n("TheRplTopologyCard"), n("TheTschScheduleCard")],
                  1
                ),
              ],
              1
            ),
            n("TheCrashReportDialog"),
          ],
          1
        );
      },
      re = [],
      le = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n("BaseKPICard", [
          n("span", { attrs: { slot: "title" }, slot: "title" }, [
            t._v("Overall E2E PDR"),
          ]),
          n(
            "span",
            { attrs: { slot: "main-contents" }, slot: "main-contents" },
            [t._v(t._s(t.pdr))]
          ),
          n(
            "div",
            { attrs: { slot: "sub-contents" }, slot: "sub-contents" },
            [
              n(
                "v-layout",
                { attrs: { "justify-center": "" } },
                [n("v-flex", { attrs: { xs10: "" } }, [n("v-divider")], 1)],
                1
              ),
              n(
                "v-layout",
                {
                  attrs: {
                    "align-center": "",
                    "justify-center": "",
                    row: "",
                    wrap: "",
                    "text-xs-center": "",
                  },
                },
                [
                  n("v-flex", { attrs: { xs3: "" } }, [
                    t._v("TX: " + t._s(t.numTx)),
                  ]),
                  n("v-flex", { attrs: { xs3: "" } }, [
                    t._v("RX: " + t._s(t.numRx)),
                  ]),
                  n("v-flex", { attrs: { xs3: "" } }, [
                    t._v("Drop: " + t._s(t.numDrop)),
                  ]),
                ],
                1
              ),
            ],
            1
          ),
          n("span", { attrs: { slot: "help" }, slot: "help" }, [
            t._v(
              "\n    End-to-End PDR since the beginning of the simulation\n  "
            ),
          ]),
        ]);
      },
      ue = [],
      ce = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs12: "", md3: "" } },
          [
            n(
              "v-tooltip",
              { attrs: { bottom: "" } },
              [
                n(
                  "v-card",
                  {
                    attrs: { slot: "activator", height: "100%" },
                    slot: "activator",
                  },
                  [
                    n(
                      "v-card-title",
                      { staticClass: "body-2 pb-0" },
                      [t._t("title")],
                      2
                    ),
                    n(
                      "v-card-text",
                      { staticClass: "title" },
                      [
                        n(
                          "v-layout",
                          { attrs: { "justify-center": "" } },
                          [t._t("main-contents")],
                          2
                        ),
                      ],
                      1
                    ),
                    t._t("sub-contents"),
                  ],
                  2
                ),
                n("span", [t._t("help")], 2),
              ],
              1
            ),
          ],
          1
        );
      },
      de = [],
      he = {},
      me = Object(f["a"])(he, ce, de, !1, null, null, null);
    me.options.__file = "BaseKPICard.vue";
    var pe = me.exports;
    v()(me, {
      VCard: x["a"],
      VCardText: b["b"],
      VCardTitle: y["a"],
      VFlex: C["a"],
      VLayout: T["a"],
      VTooltip: w["a"],
    });
    var fe = {
        components: { BaseKPICard: pe },
        mixins: [h, nt],
        data: function () {
          return { numTx: 0, numRx: 0, numDrop: 0 };
        },
        computed: {
          pdr: function () {
            return 0 === this.numTx
              ? "N/A"
              : ((this.numRx / this.numTx) * 100).toFixed(3).toString() + "%";
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.reset();
          },
          $_simulation_lastlastAppPacketEvent: function (t) {
            null !== t &&
              ((this.numTx += 1),
              "rx" === t.type ? (this.numRx += 1) : (this.numDrop += 1));
          },
        },
        methods: {
          reset: function () {
            (this.numTx = 0), (this.numRx = 0), (this.numDrop = 0);
          },
        },
      },
      _e = fe,
      ve = Object(f["a"])(_e, le, ue, !1, null, null, null);
    ve.options.__file = "TheOverallE2ePdrCard.vue";
    var ge = ve.exports;
    v()(ve, { VDivider: ft["a"], VFlex: C["a"], VLayout: T["a"] });
    var xe = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "BaseKPICard",
          [
            n("span", { attrs: { slot: "title" }, slot: "title" }, [
              t._v("Overall E2E Latency"),
            ]),
            n(
              "v-flex",
              { attrs: { slot: "main-contents" }, slot: "main-contents" },
              [
                n(
                  "v-container",
                  { attrs: { "pa-0": "", "ma-0": "" } },
                  t._l(t.items, function (e) {
                    return n(
                      "v-layout",
                      { key: e.name, attrs: { "justify-center": "" } },
                      [
                        n("v-flex", { attrs: { xs2: "", md3: "" } }, [
                          t._v(t._s(t._f("capitalize")(e.name))),
                        ]),
                        n(
                          "v-flex",
                          { attrs: { xs3: "", md4: "", "text-xs-right": "" } },
                          [t._v(t._s(t._f("latencyString")(e.value)))]
                        ),
                      ],
                      1
                    );
                  }),
                  1
                ),
              ],
              1
            ),
            n("span", { attrs: { slot: "help" }, slot: "help" }, [
              t._v(
                "\n    Average, maximum, and minimum value of End-to-End latency samples\n    since the beginning of the simulation\n  "
              ),
            ]),
          ],
          1
        );
      },
      be = [],
      ye = {
        components: { BaseKPICard: pe },
        filters: {
          capitalize: function (t) {
            return t.charAt(0).toUpperCase() + t.slice(1);
          },
          latencyString: function (t) {
            return null === t ? "N/A" : t.toFixed(2).toString() + "s";
          },
        },
        mixins: [h, nt],
        data: function () {
          return { sum: 0, numRecords: 0, max: null, min: null };
        },
        computed: {
          avg: function () {
            return 0 === this.numRecords ? null : this.sum / this.numRecords;
          },
          items: function () {
            return [
              { name: "avg", value: this.avg },
              { name: "max", value: this.max },
              { name: "min", value: this.min },
            ];
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.reset();
          },
          $_simulation_lastlastAppPacketEvent: function (t) {
            if (null !== t && ((this.numTx += 1), "rx" === t.type)) {
              var e = this.$_simulator_runningSettings.tsch_slotDuration,
                n = (t.asn - t.packet.app.timestamp) * e;
              (this.sum += n),
                (this.numRecords += 1),
                (null === this.max || this.max < n) && (this.max = n),
                (null === this.min || n < this.min) && (this.min = n);
            }
          },
        },
        methods: {
          reset: function () {
            (this.sum = 0),
              (this.numRecords = 0),
              (this.max = null),
              (this.min = null);
          },
        },
      },
      Se = ye,
      Ce = Object(f["a"])(Se, xe, be, !1, null, null, null);
    Ce.options.__file = "TheOverallE2eLatencyCard.vue";
    var $e = Ce.exports;
    v()(Ce, { VContainer: ot["a"], VFlex: C["a"], VLayout: T["a"] });
    var Te = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "BaseKPICard",
          [
            n("span", { attrs: { slot: "title" }, slot: "title" }, [
              t._v("Slotframe Occupancy"),
            ]),
            n(
              "v-flex",
              { attrs: { slot: "main-contents" }, slot: "main-contents" },
              [
                n(
                  "v-container",
                  { attrs: { "pa-0": "", "ma-0": "" } },
                  t._l(t.items, function (e) {
                    return n(
                      "v-layout",
                      { key: e.name, attrs: { "justify-center": "" } },
                      [
                        n("v-flex", { attrs: { xs2: "", md2: "" } }, [
                          t._v(t._s(t._f("capitalize")(e.name))),
                        ]),
                        n(
                          "v-flex",
                          { attrs: { xs3: "", md6: "", "text-xs-right": "" } },
                          [t._v(t._s(t._f("numCellsString")(e.value)))]
                        ),
                      ],
                      1
                    );
                  }),
                  1
                ),
              ],
              1
            ),
            n("span", { attrs: { slot: "help" }, slot: "help" }, [
              t._v(
                "\n    Maximum / minimum number of cells in one slotframe which a mote schedules\n  "
              ),
            ]),
          ],
          1
        );
      },
      we = [],
      De = n("2909"),
      ke = {
        components: { BaseKPICard: pe },
        mixins: [h, nt],
        filters: {
          capitalize: function (t) {
            return t.charAt(0).toUpperCase() + t.slice(1);
          },
          numCellsString: function (t) {
            return null === t ? "N/A" : t.toString() + " cells";
          },
        },
        data: function () {
          return { numCellsArray: [], maxNumCells: null, minNumCells: null };
        },
        computed: {
          max: function () {
            return 0 === this.numCellsArray.length
              ? null
              : Math.max.apply(Math, Object(De["a"])(this.numCellsArray));
          },
          min: function () {
            return 0 === this.numCellsArray.length
              ? null
              : Math.min.apply(Math, Object(De["a"])(this.numCellsArray));
          },
          items: function () {
            return [
              { name: "max", value: this.max },
              { name: "min", value: this.min },
            ];
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.reset();
          },
          $_simulation_lastTschCellAllocationEvent: function (t) {
            null !== t &&
              ("add" === t.type
                ? (this.numCellsArray[t.moteId] += 1)
                : "delete" === t.type && (this.numCellsArray[t.moteId] -= 1),
              (this.numCellsArray = this.numCellsArray.slice()));
          },
        },
        methods: {
          reset: function () {
            (this.numCellsArray = Array(
              this.$_simulator_runningSettings.exec_numMotes
            ).fill(0)),
              (this.maxNumCells = null),
              (this.minNumCells = null);
          },
        },
      },
      Me = ke,
      Ee = Object(f["a"])(Me, Te, we, !1, null, null, null);
    Ee.options.__file = "TheSlotframeOccupancyCard.vue";
    var Pe = Ee.exports;
    v()(Ee, { VContainer: ot["a"], VFlex: C["a"], VLayout: T["a"] });
    var Re = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "BaseKPICard",
          [
            n("span", { attrs: { slot: "title" }, slot: "title" }, [
              t._v("Mote Stats"),
            ]),
            n(
              "v-flex",
              { attrs: { slot: "main-contents" }, slot: "main-contents" },
              [
                n(
                  "v-container",
                  { attrs: { "pa-0": "", "ma-0": "" } },
                  t._l(t.items, function (e) {
                    return n(
                      "v-layout",
                      { key: e.name, attrs: { "justify-center": "" } },
                      [
                        n("v-flex", { attrs: { xs2: "", md2: "" } }, [
                          t._v(t._s(t._f("capitalize")(e.name))),
                        ]),
                        n(
                          "v-flex",
                          { attrs: { xs3: "", md7: "", "text-xs-right": "" } },
                          [t._v(t._s(e.value))]
                        ),
                      ],
                      1
                    );
                  }),
                  1
                ),
              ],
              1
            ),
            n("span", { attrs: { slot: "help" }, slot: "help" }, [
              t._v(
                "\n    The numbers of motes in the network, of synchronized motes, and of\n    securelly joined motes.\n  "
              ),
            ]),
          ],
          1
        );
      },
      Ve = [],
      Ae = {
        components: { BaseKPICard: pe },
        mixins: [h, nt],
        filters: {
          capitalize: function (t) {
            return t.charAt(0).toUpperCase() + t.slice(1);
          },
        },
        data: function () {
          return { numMotes: 0, numSynced: 0, numJoined: 0 };
        },
        computed: {
          items: function () {
            return [
              { name: "Total", value: this.numMotes },
              { name: "Synced", value: this.numSynced },
              { name: "SecJoined", value: this.numJoined },
            ];
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            "running" === t &&
              "ready" === e &&
              (this.reset(),
              (this.numMotes = this.$_simulator_runningSettings.exec_numMotes));
          },
          $_simulation_lastTschSyncEvent: function (t) {
            "tsch.synced" === t._type
              ? (this.numSynced += 1)
              : "tsch.desynced" === t._type && (this.numSynced -= 1);
          },
          $_simulation_lastSecJoinEvent: function (t) {
            "secjoin.joined" === t._type
              ? (this.numJoined += 1)
              : "secjoin.unjoined" === t._type && (this.numJoined -= 1);
          },
        },
        methods: {
          reset: function () {
            (this.numMotes = 0), (this.numSynced = 0), (this.numJoined = 0);
          },
        },
      },
      Fe = Ae,
      Oe = Object(f["a"])(Fe, Re, Ve, !1, null, null, null);
    Oe.options.__file = "TheMoteStatsCard.vue";
    var je = Oe.exports;
    v()(Oe, { VContainer: ot["a"], VFlex: C["a"], VLayout: T["a"] });
    var Ie = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "BaseChartCard",
          { attrs: { type: "line", options: t.options, series: t.series } },
          [
            n("span", { attrs: { slot: "help" }, slot: "help" }, [
              t._v(
                "\n    E2E PDR trend in the last 10 minutes. Each sample is the E2E PDR\n    value within a given minute window.\n  "
              ),
            ]),
          ]
        );
      },
      Le = [],
      Ne = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs12: "", md4: "" } },
          [
            n(
              "v-tooltip",
              { attrs: { bottom: "" } },
              [
                n(
                  "v-card",
                  {
                    attrs: { slot: "activator", height: "100%" },
                    slot: "activator",
                  },
                  [
                    n("apexchart", {
                      attrs: {
                        type: t.type,
                        height: "200px",
                        options: t.options,
                        series: t.series,
                      },
                    }),
                  ],
                  1
                ),
                n("span", [t._t("help")], 2),
              ],
              1
            ),
          ],
          1
        );
      },
      Be = [],
      He = { props: ["type", "options", "series"] },
      ze = He,
      We = Object(f["a"])(ze, Ne, Be, !1, null, null, null);
    We.options.__file = "BaseChartCard.vue";
    var Je = We.exports;
    v()(We, { VCard: x["a"], VFlex: C["a"], VTooltip: w["a"] });
    var Ue = {
        components: { BaseChartCard: Je },
        mixins: [h, nt],
        data: function () {
          return {
            options: {
              chart: {
                id: "pdr",
                group: "sync-charts",
                toolbar: { show: !1 },
                animations: { enabled: !1 },
              },
              legend: { showForSingleSeries: !0 },
              xaxis: { labels: { show: !1 } },
              yaxis: { min: 0, max: 100, tickAmount: 5 },
            },
            series: [{ name: "E2E PDR (%)", data: [{ x: 0, y: 0 }] }],
            numTx: 0,
            numRx: 0,
          };
        },
        computed: {
          pdr: function () {
            return 0 === this.numTx ? 0 : (this.numRx / this.numTx) * 100;
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.reset();
          },
          $_simulation_elapsedMinutes: function (t) {
            null !== t &&
              (this.series[0].data.push({ x: t, y: this.pdr }),
              (this.series[0].data = this.series[0].data.slice(-10)),
              (this.numTx = 0),
              (this.numRx = 0));
          },
          $_simulation_lastlastAppPacketEvent: function (t) {
            null !== t &&
              ((this.numTx += 1),
              "rx" === t.type ? (this.numRx += 1) : (this.numDrop += 1));
          },
        },
        methods: {
          reset: function () {
            (this.series[0].data = [{ x: 0, y: 0 }]),
              (this.numTx = 0),
              (this.numRx = 0);
          },
        },
      },
      Ke = Ue,
      Xe = Object(f["a"])(Ke, Ie, Le, !1, null, null, null);
    Xe.options.__file = "TheE2ePdrChartCard.vue";
    var Ge = Xe.exports,
      qe = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "BaseChartCard",
          { attrs: { type: "line", options: t.options, series: t.series } },
          [
            n("span", { attrs: { slot: "help" }, slot: "help" }, [
              t._v(
                "\n    E2E PDR trend in the last 10 minutes. Each sample is the E2E PDR\n    value within a given minute window.\n  "
              ),
            ]),
          ]
        );
      },
      Ye = [],
      Ze = {
        components: { BaseChartCard: Je },
        mixins: [h, nt],
        data: function () {
          return {
            options: {
              chart: {
                id: "latency",
                group: "sync-charts",
                toolbar: { show: !1 },
                animations: { enabled: !1 },
              },
              legend: { showForSingleSeries: !0 },
              xaxis: { labels: { show: !1 } },
              yaxis: { tickAmount: 6 },
            },
            series: [
              { name: "Avg E2E Lat (s)", data: [{ x: 0, y: 0 }] },
              { name: "Max E2E Lat (s)", data: [{ x: 0, y: 0 }] },
            ],
            numRecords: 0,
            latencySum: 0,
            latencyMax: 0,
          };
        },
        computed: {
          latencyAvg: function () {
            return 0 === this.numRecords
              ? 0
              : this.latencySum / this.numRecords;
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.reset();
          },
          $_simulation_elapsedMinutes: function (t) {
            null !== t &&
              (this.series[0].data.push({ x: t, y: this.latencyAvg }),
              (this.series[0].data = this.series[0].data.slice(-10)),
              this.series[1].data.push({ x: t, y: this.latencyMax }),
              (this.series[1].data = this.series[1].data.slice(-10)),
              (this.numRecords = 0),
              (this.latencySum = 0),
              (this.latencyMax = 0));
          },
          $_simulation_lastlastAppPacketEvent: function (t) {
            if (null !== t && "rx" === t.type) {
              var e =
                (t.asn - t.packet.app.timestamp) *
                this.$_simulator_runningSettings.tsch_slotDuration;
              (this.numRecords += 1),
                (this.latencySum += e),
                (void 0 === this.latencyMax || this.latencyMax < e) &&
                  (this.latencyMax = e);
            }
          },
        },
        methods: {
          reset: function () {
            (this.series[0].data = [{ x: 0, y: 0 }]),
              (this.series[1].data = [{ x: 0, y: 0 }]),
              (this.numRecords = 0),
              (this.latencySum = 0),
              (this.latencyMax = 0);
          },
        },
      },
      Qe = Ze,
      tn = Object(f["a"])(Qe, qe, Ye, !1, null, null, null);
    tn.options.__file = "TheE2eLatencyChartCard.vue";
    var en = tn.exports,
      nn = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "BaseChartCard",
          { attrs: { type: "bar", options: t.options, series: t.series } },
          [
            n("span", { attrs: { slot: "help" }, slot: "help" }, [
              t._v(
                "\n    Numbers of packet (unicast frame) drops with their drop reasons\n  "
              ),
            ]),
          ]
        );
      },
      an = [],
      sn = "No Route",
      on = "Full Queue",
      rn = "No Cell",
      ln = "Max Retries",
      un = "Time Exceeded",
      cn = "Rank Error",
      dn = "Full.ReassB",
      hn = "Full.VRBTbl",
      mn = {
        components: { BaseChartCard: Je },
        mixins: [h, nt],
        data: function () {
          return {
            labels: [
              { reason: "no_route", categoryLabel: sn },
              { reason: "txqueue_full", categoryLabel: on },
              { reason: "no_tx_cells", categoryLabel: rn },
              { reason: "max_retries", categoryLabel: ln },
              { reason: "time_exceeded", categoryLabel: un },
              { reason: "rank_error", categoryLabel: cn },
              { reason: "reassembly_buffer_full", categoryLabel: dn },
              { reason: "vrb_table_full", categoryLabel: hn },
            ],
            options: {
              chart: {
                id: "dropReason",
                toolbar: { show: !1 },
                animations: { enabled: !1 },
              },
              plotoptions: { distributed: !0 },
              datalabels: { enabled: !1 },
              xaxis: { categories: [sn, on, rn, ln, un, cn, dn, hn] },
              yaxis: { tickAmount: 5 },
            },
            series: [{ data: [0, 0, 0, 0, 0, 0, 0, 0] }],
          };
        },
        computed: {
          latencyAvg: function () {
            return 0 === this.numRecords
              ? 0
              : this.latencySum / this.numRecords;
          },
        },
        watch: {
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.reset();
          },
          $_simulation_lastPacketDropEvent: function (t) {
            if (null !== t) {
              var e = this.series[0].data.slice(),
                n = this.labels.find(function (e) {
                  return e.reason === t.dropReason;
                }),
                a = this.options.xaxis.categories.findIndex(function (t) {
                  return t === n.categoryLabel;
                });
              (e[a] += 1), (this.series = [{ data: e }]);
            }
          },
        },
        methods: {
          reset: function () {
            var t = Array(this.options.xaxis.categories.length).fill(0);
            this.series = [{ data: t }];
          },
        },
      },
      pn = mn,
      fn = Object(f["a"])(pn, nn, an, !1, null, null, null);
    fn.options.__file = "ThePacketDropChartCard.vue";
    var _n = fn.exports,
      vn = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs12: "", md1: "", lg4: "" } },
          [
            n(
              "v-tooltip",
              { attrs: { bottom: "" } },
              [
                n(
                  "v-card",
                  {
                    attrs: { slot: "activator", height: "100%" },
                    slot: "activator",
                  },
                  [
                    n(
                      "v-layout",
                      {
                        attrs: {
                          "align-center": "",
                          "justify-center": "",
                          row: "",
                          "fill-height": "",
                        },
                      },
                      [
                        n(
                          "v-flex",
                          [
                            n(
                              "v-container",
                              [
                                n(
                                  "v-layout",
                                  [
                                    n("cytoscape", {
                                      staticStyle: {
                                        width: "100%",
                                        height: "200px",
                                      },
                                      attrs: {
                                        config: t.config,
                                        preConfig: t.preConfig,
                                      },
                                    }),
                                  ],
                                  1
                                ),
                              ],
                              1
                            ),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
                n(
                  "span",
                  [
                    t._t("help", [
                      t._v("\n        Live view of the RPL Topology\n      "),
                    ]),
                  ],
                  2
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      gn = [],
      xn = n("b17d"),
      bn = n.n(xn),
      yn = {
        mixins: [h, nt],
        data: function () {
          return {
            options: {
              name: "dagre",
              rankDir: "TB",
              spacingFactor: 1,
              fit: !1,
            },
            config: {
              elements: [],
              layout: this.options,
              style: [
                {
                  selector: "node",
                  style: { "background-color": "#666", label: "data(id)" },
                },
                {
                  selector: "edge",
                  style: {
                    width: 3,
                    "line-color": "#ccc",
                    "target-arrow-color": "#ccc",
                    "target-arrow-shape": "triangle",
                    "curve-style": "bezier",
                  },
                },
              ],
              userZoomingEnabled: !1,
              maxZoom: 1,
            },
          };
        },
        computed: {
          breakpoint: function () {
            return this.$vuetify.breakpoint;
          },
        },
        watch: {
          breakpoint: function () {
            this.refresh();
          },
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              (this.$cytoscape.instance.then(function (t) {
                t.elements("*").remove();
              }),
              this.manipulate({ operation: "add", group: "nodes", id: 0 }),
              this.refresh());
          },
          $_simulation_lastRplParentChangeEvent: function (t) {
            null !== t &&
              (null === t.oldParentId
                ? this.manipulate({
                    operation: "add",
                    group: "nodes",
                    id: t.childId,
                  })
                : this.manipulate({
                    operation: "remove",
                    group: "edges",
                    source: t.oldParentId,
                    target: t.childId,
                  }),
              null === t.newParentId
                ? this.manipulate({
                    operation: "remove",
                    group: "nodes",
                    id: t.childId,
                  })
                : this.manipulate({
                    operation: "add",
                    group: "edges",
                    source: t.newParentId,
                    target: t.childId,
                  }),
              this.refresh());
          },
        },
        methods: {
          preConfig: function (t) {
            t.use(bn.a);
          },
          manipulate: function (t) {
            this.$cytoscape.instance.then(function (e) {
              var n = { group: t.group };
              if (
                ("nodes" === n.group
                  ? (n.id = t.id)
                  : "edges" === n.group &&
                    ((n.id = t.source + "-" + t.target),
                    (n.source = t.source),
                    (n.target = t.target)),
                "add" === t.operation)
              )
                e.add({ data: n });
              else if ("remove" === t.operation) {
                var a = e.getElementById(n.id);
                e.remove(a);
              }
            });
          },
          refresh: function () {
            var t = this;
            this.$cytoscape.instance.then(function (e) {
              e.layout(t.options).run(), e.center(), e.fit();
            });
          },
        },
      },
      Sn = yn,
      Cn = Object(f["a"])(Sn, vn, gn, !1, null, null, null);
    Cn.options.__file = "TheRplTopologyCard.vue";
    var $n = Cn.exports;
    v()(Cn, {
      VCard: x["a"],
      VContainer: ot["a"],
      VFlex: C["a"],
      VLayout: T["a"],
      VTooltip: w["a"],
    });
    var Tn = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs12: "", md11: "", lg8: "" } },
          [
            n(
              "v-tooltip",
              { attrs: { bottom: "" } },
              [
                n(
                  "v-card",
                  {
                    attrs: { slot: "activator", height: "100%" },
                    slot: "activator",
                  },
                  [
                    n(
                      "v-layout",
                      { attrs: { "align-center": "", "fill-height": "" } },
                      [
                        n(
                          "v-card-text",
                          [
                            n(
                              "v-layout",
                              {
                                attrs: {
                                  "align-center": "",
                                  "justify-center": "",
                                  row: "",
                                  wrap: "",
                                  "fill-height": "",
                                },
                              },
                              [
                                n("div", {
                                  attrs: {
                                    id: "draw-shapes",
                                    hidden: !1 === t.showMatrix,
                                  },
                                }),
                              ]
                            ),
                            t.showMatrix
                              ? n(
                                  "v-layout",
                                  {
                                    attrs: {
                                      row: "",
                                      wrap: "",
                                      "text-xs-center": "",
                                    },
                                  },
                                  [
                                    n("v-flex", { attrs: { xs12: "" } }, [
                                      n(
                                        "span",
                                        [
                                          n("font", { staticClass: "title" }, [
                                            t._v("□"),
                                          ]),
                                          t._v("No Cell, "),
                                        ],
                                        1
                                      ),
                                      n(
                                        "span",
                                        [
                                          n(
                                            "font",
                                            {
                                              staticClass: "title",
                                              attrs: { color: "blue" },
                                            },
                                            [t._v("■")]
                                          ),
                                          t._v("One TX Cell, "),
                                        ],
                                        1
                                      ),
                                      n(
                                        "span",
                                        [
                                          n(
                                            "font",
                                            {
                                              staticClass: "title",
                                              attrs: { color: "orange" },
                                            },
                                            [t._v("■")]
                                          ),
                                          t._v("Two TX Cells, "),
                                        ],
                                        1
                                      ),
                                      n(
                                        "span",
                                        [
                                          n(
                                            "font",
                                            {
                                              staticClass: "title",
                                              attrs: { color: "red" },
                                            },
                                            [t._v("■")]
                                          ),
                                          t._v(
                                            "More than two Tx Cells on that spot (SlotOffset/ChannelOffset)"
                                          ),
                                        ],
                                        1
                                      ),
                                      n("br"),
                                      t._v(
                                        "\n              (SlotOffset=0, ChannelOffset=0) is on the upper left.\n              (SlotOffset=" +
                                          t._s(this.numSlots - 1) +
                                          ", ChannelOffset=" +
                                          t._s(this.numChannels - 1) +
                                          ") is on the lower right.\n            "
                                      ),
                                    ]),
                                    n(
                                      "v-flex",
                                      { attrs: { xs12: "" } },
                                      [n("v-divider")],
                                      1
                                    ),
                                    n("v-flex", { attrs: { xs2: "" } }, [
                                      t._v("TX: " + t._s(t.numTx.value)),
                                    ]),
                                    n("v-flex", { attrs: { xs3: "" } }, [
                                      t._v(
                                        "TX/SHARED: " +
                                          t._s(t.numTxShared.value)
                                      ),
                                    ]),
                                    n("v-flex", { attrs: { xs2: "" } }, [
                                      t._v("RX: " + t._s(t.numRx.value)),
                                    ]),
                                    n("v-flex", { attrs: { xs2: "" } }, [
                                      t._v("TX/RX: " + t._s(t.numTxRx.value)),
                                    ]),
                                    n("v-flex", { attrs: { xs3: "" } }, [
                                      t._v(
                                        "TX/RX/SHARED: " +
                                          t._s(t.numTxRxShared.value)
                                      ),
                                    ]),
                                  ],
                                  1
                                )
                              : n(
                                  "v-layout",
                                  { attrs: { "justify-center": "" } },
                                  [
                                    t.showProgress
                                      ? n("v-progress-circular", {
                                          attrs: {
                                            indeterminate: "",
                                            color: "primary",
                                          },
                                        })
                                      : n(
                                          "v-flex",
                                          { attrs: { "text-xs-center": "" } },
                                          [
                                            n(
                                              "v-icon",
                                              { attrs: { color: "grey" } },
                                              [t._v("error")]
                                            ),
                                            n(
                                              "font",
                                              { attrs: { color: "grey" } },
                                              [t._v(t._s(t.errorMessage))]
                                            ),
                                          ],
                                          1
                                        ),
                                  ],
                                  1
                                ),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
                n(
                  "span",
                  [
                    t._t("help", [
                      t._v(
                        "\n        Live and global view of TSCH Schedule.\n      "
                      ),
                    ]),
                  ],
                  2
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      wn = [],
      Dn = n("53ca"),
      kn = (n("6762"), n("2fdb"), n("f077")),
      Mn = {
        mixins: [u, h, nt],
        data: function () {
          var t,
            e,
            n = 144;
          return (
            "xl" == this.$vuetify.breakpoint.name
              ? ((t = !0), (e = 1100))
              : ((t = !1), (e = 739)),
            {
              two: null,
              matrixIsAvailable: !1,
              errorMessage: null,
              cells: [],
              isLargeCanvas: t,
              canvasWidth: e,
              canvasHeight: n,
              textStyles: { size: 10 },
              textWidth: 12,
              textHeight: 12,
              minCellWidth: 3,
              minCellHeight: 3,
              minLeftPadding: 10,
              minTopPadding: 10,
              numTx: { value: 0 },
              numTxShared: { value: 0 },
              numRx: { value: 0 },
              numTxRx: { value: 0 },
              numTxRxShared: { value: 0 },
            }
          );
        },
        computed: {
          breakpoint: function () {
            return this.$vuetify.breakpoint;
          },
          minBrowserWidthToDisplay: function () {
            return this.canvasWidth + 10;
          },
          showMatrix: function () {
            return !0 === this.matrixIsAvailable && null === this.errorMessage;
          },
          showProgress: function () {
            return !1 === this.matrixIsAvailable && null === this.errorMessage;
          },
          numSlots: function () {
            return null === this.$_simulator_runningSettings
              ? 1
              : this.$_simulator_runningSettings.tsch_slotframeLength;
          },
          numChannels: function () {
            return null === this.$_simulator_runningSettings
              ? 1
              : this.$_simulator_runningSettings.phy_numChans;
          },
        },
        watch: {
          breakpoint: function () {
            !0 === this.matrixIsAvailable &&
            (this.breakpoint.width < this.minBrowserWidthToDisplay ||
              (this.isLargeCanvas && "xl" != this.$vuetify.breakpoint.name))
              ? (this.errorMessage = "Browser is too small")
              : (this.errorMessage = null);
          },
          $_simulator_operationalStatus: function (t, e) {
            (("running" === t && "ready" === e) || "aborted" === t) &&
              this.clearScheduleMatrix();
          },
          $_simulator_runningSettings: function (t) {
            null !== t && this.initializeScheduleMatrix();
          },
          $_simulation_lastTschCellAllocationEvent: function (t) {
            if (null !== t && !0 === this.matrixIsAvailable) {
              var e = this.cells[t.slotOffset][t.channelOffset];
              t.cellOptions.includes("TX") &&
                ("add" === t.type
                  ? this.updateNumTxCells(e, e.numTxCells + 1)
                  : "delete" === t.type &&
                    this.updateNumTxCells(e, e.numTxCells - 1)),
                this.updateCounter(t);
            }
          },
        },
        mounted: function () {
          this.createScheduleMatrix(),
            null !== this.$_simulator_runningSettings &&
              this.initializeScheduleMatrix();
        },
        methods: {
          createScheduleMatrix: function () {
            (this.two = new kn["a"]({
              width: this.canvasWidth,
              height: this.canvasHeight,
            }).appendTo(document.getElementById("draw-shapes"))),
              (this.cells = []);
          },
          initializeScheduleMatrix: function () {
            (this.matrixIsAvailable = !1),
              (this.errorMessage = null),
              this.resetCounters(),
              null === this.two
                ? this.createScheduleMatrix()
                : this.two.clear();
            var t = Math.floor(
                (this.canvasWidth - 2 * this.minLeftPadding) / this.numSlots
              ),
              e = Math.floor(
                (this.canvasHeight - 2 * this.minTopPadding) / this.numChannels
              );
            if ((t < e ? (e = t) : (t = e), t < this.minCellWidth))
              (this.errorMessage = "Slotframe is too long to display"),
                (this.$_app_status = "ready");
            else {
              var n =
                  (this.canvasWidth - this.textWidth - t * this.numSlots) / 2,
                a =
                  (this.canvasHeight - this.textHeight - e * this.numChannels) /
                  2;
              (this.$_app_status = "busy"),
                this.generateCells(
                  {
                    cellWidth: t,
                    cellHeight: e,
                    leftPadding: n,
                    topPadding: a,
                  },
                  { slotOffset: 0, channelOffset: 0 }
                );
            }
          },
          generateCells: function (t, e) {
            for (
              var n = this,
                a = function (a) {
                  if (a - e.slotOffset === 10) {
                    var s = n;
                    return (
                      setTimeout(function () {
                        var e = { slotOffset: a, channelOffset: 0 };
                        s.generateCells(t, e);
                      }, 5),
                      { v: void 0 }
                    );
                  }
                  n.cells[a] = [];
                  for (var i = e.channelOffset; i < n.numChannels; i++) {
                    var o = a * t.cellWidth + n.textWidth + t.leftPadding,
                      r = i * t.cellHeight + t.topPadding,
                      l = n.two.makeRectangle(o, r, t.cellWidth, t.cellHeight);
                    0 === a &&
                      i % 5 === 0 &&
                      n.two.makeText(i, t.leftPadding, r, n.textStyles),
                      a % 10 === 0 &&
                        i === n.numChannels - 1 &&
                        n.two.makeText(a, o, r + n.textHeight, n.textStyles),
                      (l.numTxCells = 0),
                      n.cells[a].push(l);
                  }
                },
                s = e.slotOffset;
              s < this.numSlots;
              s++
            ) {
              var i = a(s);
              if ("object" === Object(Dn["a"])(i)) return i.v;
            }
            this.two.update(),
              (this.matrixIsAvailable = !0),
              this.breakpoint.width < this.minBrowserWidthToDisplay
                ? (this.errorMessage = "Browser is too small")
                : (this.errorMessage = null),
              (this.$_app_status = "ready");
          },
          clearScheduleMatrix: function () {
            for (var t = 0; t < this.numSlots; t++)
              if (void 0 !== this.cells[t])
                for (var e = 0; e < this.numChannels; e++)
                  if (void 0 !== this.cells[t][e]) {
                    var n = this.cells[t][e];
                    (n.numTxCells = 0), n.noFill();
                  }
            this.two.update();
          },
          updateNumTxCells: function (t, e) {
            (t.numTxCells = e),
              0 === e
                ? t.noFill()
                : (t.fill = 1 === e ? "blue" : 2 === e ? "orange" : "red"),
              this.two.update();
          },
          updateCounter: function (t) {
            var e = !1,
              n = !1,
              a = !1,
              s = null;
            t.cellOptions.forEach(function (t) {
              "TX" === t
                ? (e = !0)
                : "RX" === t
                ? (n = !0)
                : "SHARED" === t && (a = !0);
            }),
              !0 === e && !1 === n && !1 === a
                ? (s = this.numTx)
                : !0 === e && !1 === n && !0 === a
                ? (s = this.numTxShared)
                : !1 === e && !0 === n && !1 === a
                ? (s = this.numRx)
                : !0 === e && !0 === n && !1 === a
                ? (s = this.numTxRx)
                : !0 === e && !0 === n && !0 === a && (s = this.numTxRxShared),
              "add" === t.type
                ? (s.value += 1)
                : "delete" === t.type && (s.value -= 1);
          },
          resetCounters: function () {
            (this.numTx.value = 0),
              (this.numTxShared.value = 0),
              (this.numRx.value = 0),
              (this.numTxRx.value = 0),
              (this.numTxRxShared.value = 0);
          },
        },
      },
      En = Mn,
      Pn = n("490a"),
      Rn = Object(f["a"])(En, Tn, wn, !1, null, null, null);
    Rn.options.__file = "TheTschScheduleCard.vue";
    var Vn = Rn.exports;
    v()(Rn, {
      VCard: x["a"],
      VCardText: b["b"],
      VDivider: ft["a"],
      VFlex: C["a"],
      VIcon: $["a"],
      VLayout: T["a"],
      VProgressCircular: Pn["a"],
      VTooltip: w["a"],
    });
    var An = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-dialog",
          {
            attrs: { persistent: "", width: "600px" },
            model: {
              value: t.dialog,
              callback: function (e) {
                t.dialog = e;
              },
              expression: "dialog",
            },
          },
          [
            n(
              "v-card",
              [
                n(
                  "v-flex",
                  { attrs: { xs12: "" } },
                  [
                    n(
                      "v-container",
                      { attrs: { "grid-list-md": "" } },
                      [
                        n("v-layout", {
                          attrs: { row: "", wrap: "", "fill-height": "" },
                        }),
                        n("v-card-title", { staticClass: "headline" }, [
                          t._v("Simulator Crashed :-("),
                        ]),
                        n(
                          "v-card-text",
                          [
                            n("div", { staticClass: "title" }, [
                              t._v(t._s(t.message)),
                            ]),
                            n("v-textarea", {
                              attrs: { "auto-grow": "", readonly: "" },
                              model: {
                                value: t.trace,
                                callback: function (e) {
                                  t.trace = e;
                                },
                                expression: "trace",
                              },
                            }),
                            n("div", { staticClass: "body-2" }, [
                              t._v(
                                "A crash report is available at " +
                                  t._s(t.crashReportPath)
                              ),
                            ]),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
                n(
                  "v-card-actions",
                  [
                    n("v-spacer"),
                    n(
                      "v-btn",
                      {
                        on: {
                          click: function (e) {
                            return e.stopPropagation(), t.close(e);
                          },
                        },
                      },
                      [t._v("\n        OK\n      ")]
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
          ],
          1
        );
      },
      Fn = [],
      On = {
        mixins: [h],
        computed: {
          crashReport: function () {
            return null === this.$_simulator_crashReport
              ? { trace: null, message: null, crashReportPath: null }
              : this.$_simulator_crashReport;
          },
          dialog: function () {
            return null !== this.crashReport.message;
          },
          trace: function () {
            return this.crashReport.trace;
          },
          message: function () {
            return this.crashReport.message;
          },
          crashReportPath: function () {
            return this.crashReport.crash_report_path;
          },
        },
        methods: {
          close: function () {
            this.$_simulator_clearCrashReport();
          },
        },
      },
      jn = On,
      In = n("a844"),
      Ln = Object(f["a"])(jn, An, Fn, !1, null, null, null);
    Ln.options.__file = "TheCrashReportDialog.vue";
    var Nn = Ln.exports;
    v()(Ln, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardActions: b["a"],
      VCardText: b["b"],
      VCardTitle: y["a"],
      VContainer: ot["a"],
      VDialog: S["a"],
      VFlex: C["a"],
      VLayout: T["a"],
      VSpacer: P["a"],
      VTextarea: In["a"],
    });
    var Bn = {
        components: {
          TheOverallE2ePdrCard: ge,
          TheOverallE2eLatencyCard: $e,
          TheSlotframeOccupancyCard: Pe,
          TheMoteStatsCard: je,
          TheE2ePdrChartCard: Ge,
          TheE2eLatencyChartCard: en,
          ThePacketDropChartCard: _n,
          TheRplTopologyCard: $n,
          TheTschScheduleCard: Vn,
          TheCrashReportDialog: Nn,
        },
      },
      Hn = Bn,
      zn = Object(f["a"])(Hn, oe, re, !1, null, null, null);
    zn.options.__file = "Dashboard.vue";
    var Wn = zn.exports;
    v()(zn, { VContainer: ot["a"], VFlex: C["a"], VLayout: T["a"] });
    var Jn = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs8: "" } },
          [n("TheResultStorageContainer")],
          1
        );
      },
      Un = [],
      Kn = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-container",
          { attrs: { "grid-list-sm": "", fluid: "" } },
          [
            n(
              "v-layout",
              { attrs: { "align-start": "", row: "", "fill-height": "" } },
              [
                n(
                  "v-flex",
                  { attrs: { xs12: "" } },
                  [
                    n(
                      "v-container",
                      [
                        n(
                          "v-layout",
                          { attrs: { "justify-end": "" } },
                          [
                            t.totalNumResults > 0
                              ? n(
                                  "v-btn",
                                  {
                                    attrs: { color: "error" },
                                    on: {
                                      click: function (e) {
                                        t.confirmation = !0;
                                      },
                                    },
                                  },
                                  [t._v("\n            Remove All\n          ")]
                                )
                              : t._e(),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
            0 === t.totalNumResults
              ? n(
                  "v-layout",
                  { attrs: { "justify-center": "" } },
                  [
                    n(
                      "v-flex",
                      {
                        staticClass: "title text-xs-center",
                        attrs: { x12: "" },
                      },
                      [
                        n("p", [t._v("No Result is Available under:")]),
                        n("p", [t._v(t._s(t.simDataPath))]),
                      ]
                    ),
                  ],
                  1
                )
              : t._e(),
            t._l(t.items, function (e) {
              return n(
                "v-layout",
                {
                  key: e.name,
                  attrs: {
                    "align-center": "",
                    "justify-center": "",
                    row: "",
                    wrap: "",
                    "fill-height": "",
                  },
                },
                [
                  n(
                    "v-flex",
                    { attrs: { xs12: "" } },
                    [
                      n(
                        "v-card",
                        { attrs: { hover: "" } },
                        [
                          n(
                            "v-container",
                            { attrs: { "pa-3": "" } },
                            [
                              n(
                                "v-layout",
                                {
                                  attrs: {
                                    "align-center": "",
                                    "justify-center": "",
                                    row: "",
                                  },
                                },
                                [
                                  n(
                                    "v-flex",
                                    { attrs: { xs10: "", sm3: "" } },
                                    [t._v(t._s(e.name))]
                                  ),
                                  t.$vuetify.breakpoint.smAndUp
                                    ? n(
                                        "v-flex",
                                        { attrs: { sm5: "" } },
                                        [
                                          null !== e.sfClass
                                            ? n("v-chip", [
                                                t._v(t._s(e.sfClass)),
                                              ])
                                            : t._e(),
                                          null !== e.connClass
                                            ? n("v-chip", [
                                                t._v(t._s(e.connClass)),
                                              ])
                                            : t._e(),
                                          "K7" === e.connClass
                                            ? n("v-chip", [
                                                t._v(
                                                  t._s(
                                                    t._f("traceName")(
                                                      e.connTrace
                                                    )
                                                  )
                                                ),
                                              ])
                                            : t._e(),
                                          null !== e.numMotes
                                            ? n("v-chip", [
                                                t._v(t._s(e.numMotes)),
                                              ])
                                            : t._e(),
                                        ],
                                        1
                                      )
                                    : t._e(),
                                  n("v-spacer"),
                                  t.$vuetify.breakpoint.smAndUp
                                    ? n("v-flex", { attrs: { sm2: "" } }, [
                                        t._v(
                                          "\n              " +
                                            t._s(e.lastModified) +
                                            "\n            "
                                        ),
                                      ])
                                    : t._e(),
                                  n(
                                    "v-flex",
                                    { attrs: { xs2: "", sm2: "" } },
                                    [
                                      n(
                                        "v-btn",
                                        {
                                          attrs: {
                                            href: "/result/" + e.name + ".zip",
                                            download: "",
                                            icon: "",
                                          },
                                        },
                                        [n("v-icon", [t._v("save_alt")])],
                                        1
                                      ),
                                      n(
                                        "v-btn",
                                        {
                                          attrs: { icon: "" },
                                          on: {
                                            click: function (n) {
                                              t.confirmDeletion(e.name);
                                            },
                                          },
                                        },
                                        [n("v-icon", [t._v("delete")])],
                                        1
                                      ),
                                    ],
                                    1
                                  ),
                                ],
                                1
                              ),
                            ],
                            1
                          ),
                        ],
                        1
                      ),
                    ],
                    1
                  ),
                ],
                1
              );
            }),
            n(
              "v-layout",
              { attrs: { "justify-center": "" } },
              [
                n(
                  "v-flex",
                  { attrs: { xs12: "", sm6: "" } },
                  [
                    n(
                      "v-container",
                      [
                        t.totalNumPages > 0
                          ? n(
                              "v-layout",
                              { attrs: { "justify-center": "" } },
                              [
                                n("v-pagination", {
                                  attrs: { length: t.totalNumPages },
                                  model: {
                                    value: t.currentPage,
                                    callback: function (e) {
                                      t.currentPage = e;
                                    },
                                    expression: "currentPage",
                                  },
                                }),
                              ],
                              1
                            )
                          : t._e(),
                      ],
                      1
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
            n(
              "v-dialog",
              {
                attrs: { "max-width": "500px" },
                model: {
                  value: t.confirmation,
                  callback: function (e) {
                    t.confirmation = e;
                  },
                  expression: "confirmation",
                },
              },
              [
                n(
                  "v-card",
                  [
                    n("v-card-title", { staticClass: "subheading" }, [
                      null === t.deletingItemName
                        ? n("div", [
                            t._v(
                              "\n          Delete all the files under:\n        "
                            ),
                          ])
                        : n("div", [
                            t._v(
                              "\n          Delete " +
                                t._s(t.deletingItemName) +
                                " under:\n        "
                            ),
                          ]),
                      t._v("\n          " + t._s(t.simDataPath) + "?\n      "),
                    ]),
                    n(
                      "v-layout",
                      { attrs: { "justify-center": "" } },
                      [
                        n(
                          "v-card-actions",
                          [
                            n(
                              "v-btn",
                              {
                                on: {
                                  click: function (e) {
                                    t.confirmation = !1;
                                  },
                                },
                              },
                              [t._v("Cancel")]
                            ),
                            n(
                              "v-btn",
                              {
                                attrs: { color: "error" },
                                on: { click: t.executeDeletion },
                              },
                              [t._v("Yes")]
                            ),
                          ],
                          1
                        ),
                      ],
                      1
                    ),
                  ],
                  1
                ),
              ],
              1
            ),
          ],
          2
        );
      },
      Xn = [],
      Gn = {
        filters: {
          traceName: function (t) {
            var e = t.replace(/^.*[\/]/, "");
            return e.match(/.k7.gz$/) ? e.split(".").slice(0, -2).join(".") : e;
          },
        },
        data: function () {
          return {
            currentPage: 1,
            numResultsPerPage: 5,
            totalNumResults: 0,
            results: [],
            deletingItemName: null,
            confirmation: !1,
            simDataPath: null,
          };
        },
        computed: {
          totalNumPages: function () {
            return 0 === this.numResultsPerPage
              ? 0
              : Math.floor(this.totalNumResults / this.numResultsPerPage);
          },
          now: function () {
            return Math.floor(Date.now() / 1e3);
          },
          items: function () {
            return this.results.map(function (t) {
              var e = { name: t.name, lastModified: t.last_modified },
                n = [
                  { key1: "sf_class", key2: "sfClass" },
                  { key1: "conn_class", key2: "connClass" },
                  { key1: "conn_trace", key2: "connTrace" },
                  { key1: "exec_numMotes", key2: "numMotes" },
                ];
              return (
                n.forEach(function (n) {
                  null !== t.settings && void 0 !== t.settings[n.key1]
                    ? (e[n.key2] = t.settings[n.key1])
                    : (e[n.key2] = null);
                }),
                e
              );
            });
          },
        },
        watch: {
          currentPage: function () {
            this.loadPage();
          },
          confirmation: function (t) {
            !1 === t && (this.deletingItemName = null);
          },
        },
        created: function () {
          var t = this;
          this.$_eel.get_total_number_of_results()(function (e) {
            t.totalNumResults = e;
          }),
            this.$_eel.get_sim_data_path()(function (e) {
              t.simDataPath = e;
            }),
            this.loadResults(0, this.numResultsPerPage);
        },
        methods: {
          loadPage: function () {
            this.loadResults(
              (this.currentPage - 1) * this.numResultsPerPage,
              this.numResultsPerPage
            );
          },
          loadResults: function (t, e) {
            var n = this;
            this.$_eel.get_results(
              t,
              e
            )(function (t) {
              n.results = t;
            });
          },
          deleteItem: function (t) {
            var e = this;
            this.$_eel.delete_result(t)(function () {
              (e.totalNumResults -= 1), e.loadPage();
            });
          },
          deleteAllItems: function () {
            var t = this;
            this.$_eel.delete_all_results()(function () {
              (t.totalNumResults = 0), t.loadPage(), (t.confirmation = !1);
            });
          },
          confirmDeletion: function (t) {
            (this.deletingItemName = t), (this.confirmation = !0);
          },
          executeDeletion: function () {
            null === this.deletingItemName
              ? this.deleteAllItems()
              : (this.deleteItem(this.deletingItemName),
                (this.deletingItemName = null)),
              (this.confirmation = !1);
          },
        },
      },
      qn = Gn,
      Yn = n("891e"),
      Zn = Object(f["a"])(qn, Kn, Xn, !1, null, null, null);
    Zn.options.__file = "TheResultStorageContainer.vue";
    var Qn = Zn.exports;
    v()(Zn, {
      VBtn: g["a"],
      VCard: x["a"],
      VCardActions: b["a"],
      VCardTitle: y["a"],
      VChip: it["a"],
      VContainer: ot["a"],
      VDialog: S["a"],
      VFlex: C["a"],
      VIcon: $["a"],
      VLayout: T["a"],
      VPagination: Yn["a"],
      VSpacer: P["a"],
    });
    var ta = { components: { TheResultStorageContainer: Qn } },
      ea = ta,
      na = Object(f["a"])(ea, Jn, Un, !1, null, null, null);
    na.options.__file = "Results.vue";
    var aa = na.exports;
    v()(na, { VFlex: C["a"] });
    var sa = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-flex",
          { attrs: { xs12: "" } },
          [n("TheTraceFilesContainer")],
          1
        );
      },
      ia = [],
      oa = function () {
        var t = this,
          e = t.$createElement,
          n = t._self._c || e;
        return n(
          "v-container",
          { attrs: { "grid-list-sm": "", fluid: "" } },
          [
            n(
              "v-layout",
              {
                attrs: {
                  "align-start": "",
                  "justify-center": "",
                  row: "",
                  "fill-height": "",
                },
              },
              [
                n("v-data-table", {
                  attrs: {
                    headers: t.headers,
                    items: t.files,
                    "rows-per-page-items": [
                      10,
                      25,
                      {
                        text: "$vuetify.dataIterator.rowsPerPageAll",
                        value: -1,
                      },
                    ],
                  },
                  scopedSlots: t._u([
                    {
                      key: "items",
                      fn: function (e) {
                        return [
                          n("td", [t._v(t._s(e.item.file_name))]),
                          n("td", [t._v(t._s(e.item.config.location))]),
                          n("td", [t._v(t._s(e.item.config.node_count))]),
                          n("td", [t._v(t._s(e.item.config.start_date))]),
                          n("td", [t._v(t._s(e.item.config.maxDuration))]),
                          n("td", [t._v(t._s(e.item.config.channels.length))]),
                        ];
                      },
                    },
                  ]),
                }),
              ],
              1
            ),
          ],
          1
        );
      },
      ra = [],
      la = {
        mixins: [h],
        filters: {},
        data: function () {
          return {
            headers: [
              { text: "File Name", value: "file_name" },
              { text: "Location", value: "config.location" },
              { text: "# Motes", value: "config.node_count" },
              { text: "Start Date", value: "config.start_date" },
              { text: "Duration (min)", value: "config.maxDuration" },
              { text: "# CHs", value: "config.channels.length" },
            ],
          };
        },
        computed: {
          files: function () {
            return this.$_simulator_availableTraceFiles;
          },
        },
        created: function () {},
        methods: {},
      },
      ua = la,
      ca = n("8fea"),
      da = Object(f["a"])(ua, oa, ra, !1, null, null, null);
    da.options.__file = "TheTraceFilesContainer.vue";
    var ha = da.exports;
    v()(da, { VContainer: ot["a"], VDataTable: ca["a"], VLayout: T["a"] });
    var ma = { components: { TheTraceFilesContainer: ha } },
      pa = ma,
      fa = Object(f["a"])(pa, sa, ia, !1, null, null, null);
    fa.options.__file = "TraceFiles.vue";
    var _a = fa.exports;
    v()(fa, { VFlex: C["a"] }), a["a"].use(ie["a"]);
    var va = new ie["a"]({
        mode: "history",
        base: "/",
        routes: [
          { name: "Dashboard", path: "/", component: Wn },
          { name: "Results", path: "/results", component: aa },
          { name: "Trace Files", path: "/traces", component: _a },
        ],
      }),
      ga = n("5515"),
      xa = n.n(ga);
    n("f252");
    a["a"].use(xa.a);
    var ba = n("1321"),
      ya = n.n(ba);
    a["a"].use(ya.a), a["a"].component("apexchart", ya.a);
    var Sa = {
      install: function (t) {
        t.prototype.$_eel = window.eel;
      },
    };
    a["a"].use(Sa),
      (a["a"].config.productionTip = !1),
      (window.vm = new a["a"]({
        store: se,
        router: va,
        render: function (t) {
          return t(Qt);
        },
      }).$mount("#app"));
  },
});
//# sourceMappingURL=app.1eb97c2e.js.map
