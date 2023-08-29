(() => {
    var t, n = {
        5473: (t, n, e) => {
            "use strict";
            e.d(n, {
                Z: () => i
            });
            const i = {
                APP_NAME: "Foldseek",
                APP_DESCRIPTION: "Foldseek Server offers fast and sensitive protein structure alignments against large protein structure collections",
                CITATION: 'van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, and Steinegger M. <a href="https://www.nature.com/articles/s41587-023-01773-0" target="_blank" rel="noopener">Fast and accurate protein structure search with Foldseek</a>. Nature Biotechnology, 2023.',
                NAV_URL_COUNT: "3",
                NAV_TITLE_1: "GitHub",
                NAV_URL_1: "https://foldseek.com",
                NAV_TITLE_2: "Söding Lab",
                NAV_URL_2: "https://www.mpinat.mpg.de/soeding",
                NAV_TITLE_3: "Steinegger Lab",
                NAV_URL_3: "https://steineggerlab.com/",
                QUERIES_HELP: "Enter a protein structure in PDB format or upload a PDB file.",
                UPLOAD_LABEL: "Upload PDB",
                CURL_INTRO: " Use this command to get a submit a file in PDB format to the Foldseek search server. Replace the ‘PATH_TO_FILE’ string with the path to the file.",
                MODE_HELP: "<strong>3Di/AA:</strong> fast prefilter using the 3Di alphabet and alignment using the 3Di alphabet+BLOSUM62 based Smith-Waterman-Gotoh (local alignment)<br>\n<strong>TM-align:</strong> fast prefilter using the 3Di alphabet and alignment using TM-align (global-alignment)",
                MODE_COUNT: "2",
                MODE_DEFAULT_KEY: "3diaa",
                MODE_KEY_1: "3diaa",
                MODE_TITLE_1: "3Di/AA",
                MODE_KEY_2: "tmalign",
                MODE_TITLE_2: "TM-align",
                QUERY_DEFAULT: "ATOM    866  N   PHE A 111      11.187 -12.768  -6.000\nATOM    867  CA  PHE A 111      11.895 -11.516  -5.804\nATOM    868  C   PHE A 111      13.203 -11.457  -6.592\nATOM    870  CB  PHE A 111      12.169 -11.360  -4.310\nATOM    877  N   GLY A 112      13.543 -10.277  -7.094\nATOM    878  CA  GLY A 112      14.800 -10.107  -7.788\nATOM    879  C   GLY A 112      14.816  -9.982  -9.286\nATOM    881  N   TYR A 113      13.670 -10.112  -9.938\nATOM    882  CA  TYR A 113      13.648 -10.024 -11.397\nATOM    883  C   TYR A 113      12.764  -8.904 -11.929\nATOM    885  CB  TYR A 113      13.182 -11.355 -11.997\nATOM    893  N   CYS A 114      13.052  -8.468 -13.148\nATOM    894  CA  CYS A 114      12.288  -7.406 -13.778\nATOM    895  C   CYS A 114      10.881  -7.902 -14.054\nATOM    897  CB  CYS A 114      12.938  -6.973 -15.096\nATOM    899  N   GLU A 115       9.884  -7.083 -13.740\nATOM    900  CA  GLU A 115       8.508  -7.493 -13.963\nATOM    901  C   GLU A 115       8.078  -7.419 -15.428\nATOM    903  CB  GLU A 115       7.564  -6.649 -13.087\nATOM    908  N   SER A 116       8.751  -6.604 -16.236\nATOM    909  CA  SER A 116       8.399  -6.475 -17.651\nATOM    910  C   SER A 116       9.022  -7.604 -18.460\nATOM    912  CB  SER A 116       8.874  -5.128 -18.198\nATOM    914  N   CYS A 117      10.338  -7.721 -18.376\nATOM    915  CA  CYS A 117      11.043  -8.788 -19.061\nATOM    916  C   CYS A 117      11.545  -9.657 -17.913\nATOM    918  CB  CYS A 117      12.180  -8.202 -19.896\nATOM    920  N   GLY A 118      11.749 -10.943 -18.129\nATOM    921  CA  GLY A 118      12.164 -11.781 -17.008\nATOM    922  C   GLY A 118      13.517 -11.520 -16.366\nATOM    924  N   VAL A 119      14.307 -10.654 -16.991\nATOM    925  CA  VAL A 119      15.653 -10.305 -16.546\nATOM    926  C   VAL A 119      15.839 -10.128 -15.043\nATOM    928  CB  VAL A 119      16.116  -9.004 -17.259\nATOM    931  N   GLU A 120      17.018 -10.498 -14.557\nATOM    932  CA  GLU A 120      17.318 -10.353 -13.149\nATOM    933  C   GLU A 120      17.758  -8.921 -12.929\nATOM    935  CB  GLU A 120      18.457 -11.282 -12.739\nATOM    940  N   ILE A 121      17.328  -8.325 -11.826\nATOM    941  CA  ILE A 121      17.713  -6.960 -11.477\nATOM    942  C   ILE A 121      19.000  -7.099 -10.668\nATOM    944  CB  ILE A 121      16.621  -6.291 -10.625\nATOM    948  N   GLY A 122      19.945  -6.204 -10.856\nATOM    949  CA  GLY A 122      21.175  -6.377 -10.114\nATOM    950  C   GLY A 122      21.099  -6.565  -8.605\nATOM    952  N   ILE A 123      22.051  -7.298  -8.038\nATOM    953  CA  ILE A 123      22.055  -7.474  -6.607\nATOM    954  C   ILE A 123      22.389  -6.135  -5.992\nATOM    956  CB  ILE A 123      23.078  -8.512  -6.173\nATOM    960  N   ARG A 124      23.412  -5.481  -6.521\nATOM    961  CA  ARG A 124      23.804  -4.174  -5.993\nATOM    962  C   ARG A 124      22.719  -3.163  -6.291\nATOM    964  CB  ARG A 124      25.110  -3.680  -6.625\nATOM    971  N   ARG A 125      21.969  -3.378  -7.364\nATOM    972  CA  ARG A 125      20.903  -2.436  -7.674\nATOM    973  C   ARG A 125      19.754  -2.611  -6.682\nATOM    975  CB  ARG A 125      20.358  -2.617  -9.083\nATOM    982  N   LEU A 126      19.493  -3.856  -6.289\nATOM    983  CA  LEU A 126      18.430  -4.140  -5.333\nATOM    984  C   LEU A 126      18.838  -3.655  -3.951\nATOM    986  CB  LEU A 126      18.141  -5.637  -5.271\nATOM    990  N   GLU A 127      20.138  -3.596  -3.708\nATOM    991  CA  GLU A 127      20.632  -3.131  -2.429\nATOM    992  C   GLU A 127      20.396  -1.621  -2.356\nATOM    994  CB  GLU A 127      22.117  -3.451  -2.320\nATOM    999  N   ALA A 128      20.326  -0.979  -3.520\nATOM   1000  CA  ALA A 128      20.074   0.459  -3.603\nATOM   1001  C   ALA A 128      18.574   0.724  -3.409\nATOM   1003  CB  ALA A 128      20.517   0.985  -4.943\nATOM   1004  N   ARG A 129      17.730   0.026  -4.174\nATOM   1005  CA  ARG A 129      16.277   0.152  -4.044\nATOM   1006  C   ARG A 129      15.726  -1.263  -4.110\nATOM   1008  CB  ARG A 129      15.680   0.998  -5.173\nATOM   1015  N   PRO A 130      15.684  -1.961  -2.968\nATOM   1016  CA  PRO A 130      15.183  -3.334  -2.892\nATOM   1017  C   PRO A 130      13.742  -3.504  -3.336\nATOM   1019  CB  PRO A 130      15.393  -3.691  -1.429\nATOM   1022  N   THR A 131      13.075  -2.383  -3.540\nATOM   1023  CA  THR A 131      11.675  -2.355  -3.940\nATOM   1024  C   THR A 131      11.531  -2.277  -5.471\nATOM   1026  CB  THR A 131      11.004  -1.137  -3.239\nATOM   1029  N   ALA A 132      12.661  -2.293  -6.172\nATOM   1030  CA  ALA A 132      12.672  -2.208  -7.625\nATOM   1031  C   ALA A 132      11.798  -3.246  -8.352\nATOM   1033  CB  ALA A 132      14.106  -2.304  -8.114\nATOM   1034  N   ASP A 133      10.971  -2.777  -9.287\nATOM   1035  CA  ASP A 133      10.071  -3.635 -10.060\nATOM   1036  C   ASP A 133      10.581  -3.912 -11.473\nATOM   1038  CB  ASP A 133       8.681  -2.987 -10.220\nATOM   1042  N   LEU A 134      11.366  -2.982 -12.010\nATOM   1043  CA  LEU A 134      11.863  -3.127 -13.369\nATOM   1044  C   LEU A 134      13.361  -3.082 -13.523\nATOM   1046  CB  LEU A 134      11.257  -2.039 -14.242\nATOM   1050  N   CYS A 135      13.836  -3.733 -14.589\nATOM   1051  CA  CYS A 135      15.243  -3.648 -14.882\nATOM   1052  C   CYS A 135      15.282  -2.173 -15.324\nATOM   1054  CB  CYS A 135      15.651  -4.622 -16.008\nATOM   1056  N   ILE A 136      16.461  -1.566 -15.338\nATOM   1057  CA  ILE A 136      16.567  -0.158 -15.714\nATOM   1058  C   ILE A 136      15.950   0.181 -17.061\nATOM   1060  CB  ILE A 136      18.043   0.319 -15.697\nATOM   1064  N   ASP A 137      16.145  -0.690 -18.047\nATOM   1065  CA  ASP A 137      15.602  -0.448 -19.378\nATOM   1066  C   ASP A 137      14.082  -0.391 -19.394\nATOM   1068  CB  ASP A 137      16.048  -1.516 -20.372\nATOM   1072  N   CYS A 138      13.433  -1.411 -18.854\nATOM   1073  CA  CYS A 138      11.977  -1.428 -18.842\nATOM   1074  C   CYS A 138      11.458  -0.325 -17.968\nATOM   1076  CB  CYS A 138      11.431  -2.759 -18.330\nATOM   1078  N   LYS A 139      12.159  -0.068 -16.872\nATOM   1079  CA  LYS A 139      11.752   0.988 -15.957\nATOM   1080  C   LYS A 139      11.752   2.318 -16.682\nATOM   1082  CB  LYS A 139      12.709   1.093 -14.766\nATOM   1087  N   THR A 140      12.841   2.584 -17.394\nATOM   1088  CA  THR A 140      12.987   3.830 -18.134\nATOM   1089  C   THR A 140      12.001   3.945 -19.284\nATOM   1091  CB  THR A 140      14.413   3.980 -18.671\nATOM   1094  N   LEU A 141      11.855   2.866 -20.038\nATOM   1095  CA  LEU A 141      10.936   2.857 -21.156\nATOM   1096  C   LEU A 141       9.543   3.165 -20.663\nATOM   1098  CB  LEU A 141      10.967   1.509 -21.855\nATOM   1102  N   ALA A 142       9.202   2.630 -19.501\nATOM   1103  CA  ALA A 142       7.888   2.875 -18.910\nATOM   1104  C   ALA A 142       7.720   4.354 -18.613\nATOM   1106  CB  ALA A 142       7.734   2.069 -17.624\nATOM   1107  N   GLU A 143       8.760   4.969 -18.070\nATOM   1108  CA  GLU A 143       8.715   6.382 -17.737\nATOM   1109  C   GLU A 143       8.556   7.223 -18.995\nATOM   1111  CB  GLU A 143       9.992   6.783 -17.003\nATOM   1116  N   ILE A 144       9.188   6.790 -20.080\nATOM   1117  CA  ILE A 144       9.096   7.513 -21.329\nATOM   1118  C   ILE A 144       7.684   7.397 -21.873\nATOM   1120  CB  ILE A 144      10.091   6.976 -22.380\nATOM   1124  N   ARG A 145       7.153   6.178 -21.916\nATOM   1125  CA  ARG A 145       5.798   5.945 -22.417\nATOM   1126  C   ARG A 145       4.846   6.844 -21.651\nATOM   1128  CB  ARG A 145       5.359   4.495 -22.200\nATOM   1135  N   GLU A 146       5.063   6.922 -20.346\nATOM   1136  CA  GLU A 146       4.263   7.735 -19.443\nATOM   1137  C   GLU A 146       4.121   9.167 -19.951\nATOM   1139  CB  GLU A 146       4.936   7.716 -18.080\nATOM   1144  N   LYS A 147       5.248   9.860 -20.097\nATOM   1145  CA  LYS A 147       5.253  11.240 -20.581\nATOM   1146  C   LYS A 147       4.540  11.421 -21.924\nATOM   1148  CB  LYS A 147       6.693  11.757 -20.710\nATOM   1153  N   GLN A 148       4.576  10.393 -22.762\nATOM   1154  CA  GLN A 148       3.951  10.453 -24.085\nATOM   1155  C   GLN A 148       2.471  10.044 -24.106\nATOM   1157  CB  GLN A 148       4.750   9.592 -25.070\nATOM   1162  N   MET A 149       2.128   8.997 -23.359\nATOM   1163  CA  MET A 149       0.743   8.529 -23.282\nATOM   1164  C   MET A 149      -0.049   9.525 -22.433\nATOM   1166  CB  MET A 149       0.660   7.141 -22.624\nATOM   1170  N   ALA A 150       0.664  10.279 -21.603\nATOM   1171  CA  ALA A 150       0.044  11.272 -20.740\nATOM   1172  C   ALA A 150      -0.134  12.585 -21.497\nATOM   1174  CB  ALA A 150       0.902  11.499 -19.503\nATOM   1175  N   GLY A 151       0.960  13.327 -21.647\nATOM   1176  CA  GLY A 151       0.909  14.596 -22.353\nATOM   1177  C   GLY A 151       0.566  14.495 -23.835\nTER"
            };
        },
        8615: (t, n, e) => {
            "use strict";
            e.d(n, {
                Z: () => i
            });
            const i = {
                APP_NAME: "MMseqs2",
                APP_DESCRIPTION: "MMseqs2 server offers fast and sensitive sequence alignments against large sequence databases",
                CITATION: "Mirdita M., Steinegger M., and Söding J., <a href=“https://doi.org/10.1093/bioinformatics/bty1057” target=“_blank” rel=“noopener”>MMseqs2 desktop and local web server app for fast, interactive sequence searches</a>, <i>Bioinformatics</i>, 2019.",
                NAV_URL_COUNT: "3",
                NAV_TITLE_1: "MMseqs2",
                NAV_URL_1: "https://mmseqs.com",
                NAV_TITLE_2: "GitHub",
                NAV_URL_2: "https://github.com/soedinglab/MMseqs2-App",
                NAV_TITLE_3: "Södinglab",
                NAV_URL_3: "http://www.mpibpc.mpg.de/soeding",
                QUERIES_HELP: "Enter a list of either protein or nucleotide sequences in FASTA format or upload a FASTA file.",
                UPLOAD_LABEL: "Upload FASTA File",
                CURL_INTRO: " Use this command to get a submit a file in fasta format to the MMseqs2 search server. Replace the ‘PATH_TO_FILE’ string with the path to the file.",
                MODE_HELP: "‘All’ shows all hits under an e-value cutoff. ‘Greedy Best Hits’ tries to cover the search query.",
                MODE_COUNT: "2",
                MODE_DEFAULT_KEY: "accept",
                MODE_KEY_1: "accept",
                MODE_TITLE_1: "All Hits",
                MODE_KEY_2: "summary",
                MODE_TITLE_2: "Greedy Best Hits",
                MODE_KEY_3: "",
                MODE_TITLE_3: "",
                QUERY_DEFAULT: ">TEST\nMPKIIEAIYENGVFKPLQKVDLKEGEKAKIVLESISDKTFGILKASETEIKKVLEEIDDFWGVC"
            };
        },
        8217: (t, n, e) => {
            "use strict";
            var i = e(144), a = e(6828), r = e(1002), s = {
                selector: "vue-portal-target"
            };
            const o = s;
            var l = "undefined" != typeof window && void 0 !== ("undefined" == typeof document ? "undefined" : (0, 
            r.Z)(document));
            const A = i.Z.extend({
                abstract: !0,
                name: "PortalOutlet",
                props: [ "nodes", "tag" ],
                data: function(t) {
                    return {
                        updatedNodes: t.nodes
                    };
                },
                render: function(t) {
                    var n = this.updatedNodes && this.updatedNodes();
                    return n ? 1 !== n.length || n[0].text ? t(this.tag || "DIV", n) : n : t();
                },
                destroyed: function() {
                    var t = this.$el;
                    t && t.parentNode.removeChild(t);
                }
            }), c = i.Z.extend({
                name: "VueSimplePortal",
                props: {
                    disabled: {
                        type: Boolean
                    },
                    prepend: {
                        type: Boolean
                    },
                    selector: {
                        type: String,
                        default: function() {
                            return "#".concat(o.selector);
                        }
                    },
                    tag: {
                        type: String,
                        default: "DIV"
                    }
                },
                render: function(t) {
                    if (this.disabled) {
                        var n = this.$scopedSlots && this.$scopedSlots.default();
                        return n ? n.length < 2 && !n[0].text ? n : t(this.tag, n) : t();
                    }
                    return t();
                },
                created: function() {
                    this.getTargetEl() || this.insertTargetEl();
                },
                updated: function() {
                    var t = this;
                    this.$nextTick((function() {
                        t.disabled || t.slotFn === t.$scopedSlots.default || (t.container.updatedNodes = t.$scopedSlots.default), 
                        t.slotFn = t.$scopedSlots.default;
                    }));
                },
                beforeDestroy: function() {
                    this.unmount();
                },
                watch: {
                    disabled: {
                        immediate: !0,
                        handler: function(t) {
                            t ? this.unmount() : this.$nextTick(this.mount);
                        }
                    }
                },
                methods: {
                    getTargetEl: function() {
                        if (l) return document.querySelector(this.selector);
                    },
                    insertTargetEl: function() {
                        if (l) {
                            var t = document.querySelector("body"), n = document.createElement(this.tag);
                            n.id = this.selector.substring(1), t.appendChild(n);
                        }
                    },
                    mount: function() {
                        if (l) {
                            var t = this.getTargetEl(), n = document.createElement("DIV");
                            this.prepend && t.firstChild ? t.insertBefore(n, t.firstChild) : t.appendChild(n), 
                            this.container = new A({
                                el: n,
                                parent: this,
                                propsData: {
                                    tag: this.tag,
                                    nodes: this.$scopedSlots.default
                                }
                            });
                        }
                    },
                    unmount: function() {
                        this.container && (this.container.$destroy(), delete this.container);
                    }
                }
            });
            function d(t) {
                var n, e = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : {};
                t.component(e.name || "portal", c), e.defaultSelector && (n = e.defaultSelector, 
                s.selector = n);
            }
            "undefined" != typeof window && window.Vue && window.Vue === i.Z && i.Z.use(d);
            const u = d;
            var h = e(5317), p = (e(8197), e(7895), e(1434), function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("v-app", {
                    class: {
                        electron: t.$ELECTRON
                    },
                    attrs: {
                        id: "app"
                    }
                }, [ e("v-main", [ e("ResultLocal") ], 1) ], 1);
            });
            p._withStripped = !0;
            var g = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("div", [ e("v-app-bar", {
                    attrs: {
                        app: "",
                        height: "48px",
                        fixed: "",
                        "clipped-left": ""
                    }
                }, [ e("img", {
                    attrs: {
                        height: "28px",
                        src: "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbDpzcGFjZT0icHJlc2VydmUiIHN0eWxlPSJmaWxsLXJ1bGU6ZXZlbm9kZDtjbGlwLXJ1bGU6ZXZlbm9kZDtzdHJva2UtbGluZWNhcDpyb3VuZDtzdHJva2UtbGluZWpvaW46cm91bmQ7c3Ryb2tlLW1pdGVybGltaXQ6MTAiIHZpZXdCb3g9IjAgMCA0NjggMzA2Ij48cGF0aCBkPSJNMzcyIDIwMnMxNC0xIDM3LTE5YzIzLTE3IDQwLTQ5IDU1LTU1bC0xMTQgMjQtNCAzMiAyNiAxOFoiIHN0eWxlPSJmaWxsOiNmN2QxOGE7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik02MiAxMzlTODcgMjEgMjY5IDJsMSAxLTQ2IDYxcy00MC0zLTU1IDdjMCAwIDE5LTEzIDY5LTRzNTAtMjAgNTAtMjAgOCAyMiAwIDI5bDI5IDE0LTE4IDRzMTI1LTEyIDE2NyAzM2MwIDAtMjYgMTctNjAgMjAtNTYgNS02MiAyMi02MiAyMnMyNS0xMCA0MyA0bC0yMiA5czE1IDggMTUgMjNsLTI2IDEwczM2LTE4IDUyLTdsLTI0IDE4czIzIDMgMzggMTVsLTMyIDhzMTUgMiAyNyAzMWwtNDUtNnM3IDkgNCAzMGwtMjUtMjJzLTE3IDQ2LTE1OCAyQzQ5IDI0MCA1NiAyMjEgNTAgMTkxbC0yNi0xczItMTUgMTgtMjFMMiAxNDJzMjQtMTMgNDItOGwtOC0yNXMyOSAxMSAyNiAzMFoiIHN0eWxlPSJmaWxsOiNlMTMyMTM7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0xMDEgMjUzYy00Ni0yMyA4LTEzNCAzNy0xNTEgMjgtMTYgNTcgNyA2MyAxOSAwIDAgMjMtMTggNTctN3M0OSA0NyAzNiAxMTVjLTggNDEtMjQgNTgtMzUgNjUtNyA0LTE0IDUtMjEgMy0yNS02LTEwNS0yNy0xMzctNDRaIiBzdHlsZT0iZmlsbDojZjdkMThhO2ZpbGwtcnVsZTpub256ZXJvO3N0cm9rZTojMDAwO3N0cm9rZS13aWR0aDo0LjQ4cHgiLz48cGF0aCBkPSJNMTM2IDExMnMtNDEtMTAtNTYgMThjLTE1IDI3IDEyIDM4IDI3IDQzIDE2IDQgNDcgNCA1Ny0xM3MtMS0zOC0yOC00OFoiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0xMTYgMTYwYzE2IDggMzQtMzcgMjAtNDQtMTQtNi00MCAzNS0yMCA0NFoiIHN0eWxlPSJmaWxsLXJ1bGU6bm9uemVybztzdHJva2U6IzAwMDtzdHJva2Utd2lkdGg6NC40OHB4Ii8+PHBhdGggZD0iTTI4NCAxNDhjLTQxLTE1LTU5IDUtNjUgMjJzMiA0NCA0MiA1MyA1MC00IDU2LTE5YzUtMTYgNi00MS0zMy01NloiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0yNDggMTk5YzE5IDkgNDctNDEgMjMtNTJzLTQzIDQzLTIzIDUyWm0tODUtMTVjMS04IDIwLTEgMjAgNSAwIDctOSA4LTEyIDctNC0xLTktNi04LTEyWiIgc3R5bGU9ImZpbGwtcnVsZTpub256ZXJvO3N0cm9rZTojMDAwO3N0cm9rZS13aWR0aDo0LjQ4cHgiLz48cGF0aCBkPSJNMTMyIDEyMGM3IDMtMiAxNS02IDEyczMtMTQgNi0xMlptMTI4IDMwYzcgMy0yIDE1LTYgMTItNC0yIDMtMTQgNi0xMloiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm8iLz48cGF0aCBkPSJtMTE1IDIxMiA5LTRzLTggNyAwIDEzYzggNyAyNS00IDQ2LTEgMjEgNCA0MCAxOSA1NSAyMSAxNiAzIDI0IDEgMjMtNC0xLTYgNSA3IDUgNyIgc3R5bGU9ImZpbGw6bm9uZTtmaWxsLXJ1bGU6bm9uemVybztzdHJva2U6IzAwMDtzdHJva2Utd2lkdGg6NC40OHB4Ii8+PC9zdmc+"
                    }
                }), t._v("\n         \n        "), e("v-app-bar-title", {
                    staticClass: "ml-2"
                }, [ t._v(t._s(t.$STRINGS.APP_NAME) + " Search") ]), t._v(" "), e("v-spacer"), t._v(" "), e("v-file-input", {
                    staticClass: "shrink",
                    staticStyle: {
                        position: "relative",
                        top: "30%"
                    },
                    attrs: {
                        id: "uploadData",
                        type: "file",
                        accept: "application/json",
                        placeholder: "Load JSON data file",
                        "single-line": "",
                        outlined: "",
                        filled: "",
                        flat: "",
                        dense: ""
                    },
                    on: {
                        change: t.uploadData
                    }
                }), t._v(" "), e("v-toolbar-items", [ e("v-btn", {
                    attrs: {
                        text: ""
                    },
                    on: {
                        click: t.downloadData
                    }
                }, [ e("v-icon", [ t._v("\n                    " + t._s(t.$MDI.FileDownloadOutline) + "\n                ") ]) ], 1), t._v(" "), t._l(t.$STRINGS.NAV_URL_COUNT - 0, (function(n) {
                    return e("v-btn", {
                        key: n,
                        staticClass: "hidden-sm-and-down",
                        attrs: {
                            text: "",
                            rel: "external noopener",
                            target: "_blank",
                            href: t.$STRINGS["NAV_URL_" + n]
                        }
                    }, [ t._v(t._s(t.$STRINGS["NAV_TITLE_" + n])) ]);
                })) ], 2) ], 1), t._v(" "), t.hits ? e("v-tabs", {
                    staticStyle: {
                        "margin-bottom": "1em"
                    },
                    attrs: {
                        "center-active": "",
                        grow: "",
                        "show-arrows": ""
                    }
                }, t._l(t.hits, (function(n, i) {
                    return e("v-tab", {
                        key: n.query.header,
                        on: {
                            click: function(n) {
                                return t.changeResult(i);
                            }
                        }
                    }, [ t._v("\n            " + t._s(n.query.header) + " (" + t._s(n.results[0].alignments ? n.results[0].alignments.length : 0) + ")\n        ") ]);
                })), 1) : t._e(), t._v(" "), t.hits ? e("ResultView", {
                    key: t.currentIndex,
                    attrs: {
                        ticket: t.ticket,
                        error: t.error,
                        mode: t.mode,
                        hits: t.currentResult,
                        selectedDatabases: t.selectedDatabases,
                        tableMode: t.tableMode
                    }
                }) : e("v-container", {
                    attrs: {
                        "grid-list-md": "",
                        fluid: "",
                        "pa-2": ""
                    }
                }, [ e("v-layout", {
                    attrs: {
                        wrap: ""
                    }
                }, [ e("v-flex", {
                    attrs: {
                        xs12: ""
                    }
                }, [ e("v-card", {
                    attrs: {
                        rounded: "0"
                    }
                }, [ e("v-card-title", {
                    staticClass: "mb-0 pa-4",
                    attrs: {
                        "primary-title": ""
                    }
                }, [ t._v("\n                        No data loaded\n                    ") ]) ], 1) ], 1) ], 1) ], 1), t._v(" "), e("v-container", {
                    attrs: {
                        "grid-list-md": "",
                        fluid: "",
                        "pa-2": ""
                    }
                }, [ e("v-layout", {
                    attrs: {
                        wrap: ""
                    }
                }, [ e("v-flex", {
                    attrs: {
                        xs12: ""
                    }
                }, [ e("v-card", {
                    attrs: {
                        rounded: "0"
                    }
                }, [ e("v-card-title", {
                    staticClass: "pb-0 mb-0",
                    attrs: {
                        "primary-title": ""
                    }
                }, [ e("div", {
                    staticClass: "text-h5 mb-0"
                }, [ t._v("Reference") ]) ]), t._v(" "), e("v-card-title", {
                    staticClass: "pt-0 mt-0",
                    attrs: {
                        "primary-title": ""
                    }
                }, [ e("p", {
                    staticClass: "text-subtitle-2 mb-0",
                    domProps: {
                        innerHTML: t._s(t.$STRINGS.CITATION)
                    }
                }) ]) ], 1) ], 1) ], 1) ], 1) ], 1);
            };
            g._withStripped = !0;
            var m = "1f77b4aec7e8ff7f0effbb782ca02c98df8ad62728ff98969467bdc5b0d58c564bc49c94e377c2f7b6d27f7f7fc7c7c7bcbd22dbdb8d17becf9edae5".match(/.{6}/g).map((function(t) {
                return "#" + t;
            }));
            function v(t) {
                t = function(t) {
                    var n = function(t) {
                        return parseInt(t, 16) / 255;
                    };
                    return [ n(t.slice(1, 3)), n(t.slice(3, 5)), n(t.slice(5, 7)) ];
                }(t);
                var n = t[0], e = t[1], i = t[2], a = Math.min(n, e, i), r = Math.max(n, e, i), s = NaN, o = r - a, l = (r + a) / 2;
                return o ? (s = n === r ? (e - i) / o + 6 * (e < i) : e === r ? (i - n) / o + 2 : (n - e) / o + 4, 
                o /= l < .5 ? r + a : 2 - r - a, s *= 60) : o = l > 0 && l < 1 ? 0 : s, [ s, o, l ];
            }
            function f(t, n) {
                var e = "undefined" != typeof Symbol && t[Symbol.iterator] || t["@@iterator"];
                if (!e) {
                    if (Array.isArray(t) || (e = function(t, n) {
                        if (!t) return;
                        if ("string" == typeof t) return b(t, n);
                        var e = Object.prototype.toString.call(t).slice(8, -1);
                        "Object" === e && t.constructor && (e = t.constructor.name);
                        if ("Map" === e || "Set" === e) return Array.from(t);
                        if ("Arguments" === e || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(e)) return b(t, n);
                    }(t)) || n && t && "number" == typeof t.length) {
                        e && (t = e);
                        var i = 0, a = function() {};
                        return {
                            s: a,
                            n: function() {
                                return i >= t.length ? {
                                    done: !0
                                } : {
                                    done: !1,
                                    value: t[i++]
                                };
                            },
                            e: function(t) {
                                throw t;
                            },
                            f: a
                        };
                    }
                    throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.");
                }
                var r, s = !0, o = !1;
                return {
                    s: function() {
                        e = e.call(t);
                    },
                    n: function() {
                        var t = e.next();
                        return s = t.done, t;
                    },
                    e: function(t) {
                        o = !0, r = t;
                    },
                    f: function() {
                        try {
                            s || null == e.return || e.return();
                        } finally {
                            if (o) throw r;
                        }
                    }
                };
            }
            function b(t, n) {
                (null == n || n > t.length) && (n = t.length);
                for (var e = 0, i = new Array(n); e < n; e++) i[e] = t[e];
                return i;
            }
            const C = {
                name: "result",
                data: function() {
                    return {
                        ticket: "",
                        error: "",
                        mode: "",
                        hits: null,
                        alignment: null,
                        activeTarget: null,
                        alnBoxOffset: 0,
                        selectedDatabases: 0,
                        tableMode: 0
                    };
                },
                methods: {
                    resetProperties: function() {},
                    fetchData: function() {},
                    tryLinkTargetToDB: function(t, n) {
                        var e = n.toLowerCase();
                        return e.startsWith("pfam") ? "https://pfam.xfam.org/family/" + t : e.startsWith("pdb") ? "https://www.rcsb.org/pdb/explore.do?structureId=" + t.replaceAll(/\.(cif|pdb)(\.gz)?/g, "").split("_")[0] : e.startsWith("uniclust") || e.startsWith("uniprot") || e.startsWith("sprot") || e.startsWith("swissprot") ? "https://www.uniprot.org/uniprot/" + t : e.startsWith("eggnog_") ? "http://eggnogdb.embl.de/#/app/results?target_nogs=" + t : e.startsWith("cdd") ? "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=" + t : t.startsWith("AF-") ? "https://www.alphafold.ebi.ac.uk/entry/" + t.replaceAll(/-F[0-9]+-model_v[0-9]+(\.(cif|pdb))?(\.gz)?(_[A-Z0-9]+)?$/g, "") : t.startsWith("GMGC") ? "https://gmgc.embl.de/search.cgi?search_id=" + t.replaceAll(/\.(cif|pdb)(\.gz)?/g, "") : t.startsWith("MGYP") ? "https://esmatlas.com/explore/detail/" + t.replaceAll(/\.(cif|pdb)(\.gz)?/g, "") : e.startsWith("cath") ? t.startsWith("af_") ? "https://www.cathdb.info/version/latest/superfamily/" + t.substring(t.lastIndexOf("_") + 1) : "https://www.cathdb.info/version/latest/domain/" + t : null;
                    },
                    tryFixTargetName: function(t, n) {
                        var e = n.toLowerCase();
                        if (t.startsWith("AF-")) return t.replaceAll(/\.(cif|pdb)(\.gz)?(_[A-Z0-9]+)?$/g, "");
                        if (e.startsWith("pdb") || e.startsWith("gmgc") || e.startsWith("mgyp") || e.startsWith("mgnify")) return t.replaceAll(/\.(cif|pdb)(\.gz)?/g, "");
                        if (e.startsWith("cath") && t.startsWith("af_")) {
                            var i = t.match(/^af_([A-Z0-9]+)_(\d+)_(\d+)_(\d+\.\d+\.\d+\.\d+)$/);
                            if (i && 5 == i.length) return i[4] + " " + i[1] + " " + i[2] + "-" + i[3];
                        }
                        return t;
                    },
                    setColorScheme: function() {
                        if (this.hits) {
                            var t, n, e, i, a, r, s, o, l, A = (t = [], n = 1, function(e) {
                                var i = e + "", a = t[i];
                                return a || (a = t[i] = n++), m[(a - 1) % m.length];
                            }), c = f(this.currentResult.results);
                            try {
                                for (c.s(); !(e = c.n()).done; ) {
                                    var d = e.value;
                                    d.color = A(d.db);
                                    var u, h = v(d.color), p = {
                                        score: Number.MIN_VALUE
                                    }, g = {
                                        score: Number.MAX_VALUE
                                    }, b = f(d.alignments);
                                    try {
                                        for (b.s(); !(u = b.n()).done; ) {
                                            var C = u.value;
                                            for (var M in g) g[M] = C[M] < g[M] ? C[M] : g[M], p[M] = C[M] > p[M] ? C[M] : p[M];
                                        }
                                    } catch (t) {
                                        b.e(t);
                                    } finally {
                                        b.f();
                                    }
                                    var y, w = f(d.alignments);
                                    try {
                                        for (w.s(); !(y = w.n()).done; ) {
                                            var x = y.value, I = (s = g.score / p.score, o = 1, l = x.score / p.score, s * (1 - l) + o * l), S = (i = h[2] * Math.pow(.55, -(1 - I)), 
                                            a = .1, r = .9, Math.max(a, Math.min(r, i)));
                                            x.color = "hsl(".concat(h[0], ", ").concat(100 * h[1], "%, ").concat(100 * S, "%)");
                                        }
                                    } catch (t) {
                                        w.e(t);
                                    } finally {
                                        w.f();
                                    }
                                }
                            } catch (t) {
                                c.e(t);
                            } finally {
                                c.f();
                            }
                        }
                    },
                    parseResults: function(t) {
                        var n = 0, e = 0;
                        for (var i in t.results) {
                            var a = t.results[i], r = a.db;
                            for (var s in a.hasDescription = !1, a.hasTaxonomy = !1, null == a.alignments && n++, 
                            e++, a.alignments) {
                                var o = a.alignments[s], l = o.target.split(" ");
                                o.target = l[0], o.description = l.slice(1).join(" "), o.description.length > 1 && (a.hasDescription = !0), 
                                o.href = this.tryLinkTargetToDB(o.target, r), o.target = this.tryFixTargetName(o.target, r), 
                                o.id = "result-" + i + "-" + s, o.active = !1, "tmalign" != this.mode && (o.eval = "string" == typeof o.eval ? o.eval : o.eval.toExponential(2)), 
                                o.prob = "string" == typeof o.prob ? o.prob : o.prob.toFixed(2), "tmalign" == this.mode && (o.eval = "string" == typeof o.eval ? o.eval : o.eval.toFixed(3)), 
                                "taxId" in o && (a.hasTaxonomy = !0);
                            }
                        }
                        return 0 != e && n / e == 1 ? {
                            results: []
                        } : t;
                    }
                },
                watch: {
                    hits: function() {
                        this.setColorScheme();
                    }
                }
            };
            var M = e(1900), y = (0, M.Z)(C, undefined, undefined, !1, null, null, null);
            y.options.__file = "frontend/ResultMixin.vue";
            const w = y.exports;
            var x = function() {
                var t = this, n = t.$createElement, i = t._self._c || n;
                return i("v-container", {
                    attrs: {
                        "grid-list-md": "",
                        fluid: "",
                        "pa-2": ""
                    }
                }, [ i("v-layout", {
                    attrs: {
                        wrap: ""
                    }
                }, [ i("v-flex", {
                    attrs: {
                        xs12: ""
                    }
                }, [ i("panel", [ i("template", {
                    slot: "header"
                }, [ t.$LOCAL || t.hits && t.hits.query ? [ i("span", {
                    staticClass: "hidden-sm-and-down"
                }, [ t._v("Results: ") ]), t._v(" "), i("small", {
                    staticClass: "ticket"
                }, [ t._v(t._s(t.hits.query.header)) ]) ] : [ i("span", {
                    staticClass: "hidden-sm-and-down"
                }, [ t._v("Results for job: ") ]), t._v(" "), i("small", {
                    staticClass: "ticket"
                }, [ t._v(t._s(t.ticket)) ]) ] ], 2), t._v(" "), t.$LOCAL || "PENDING" != t.resultState ? t.$LOCAL || "EMPTY" != t.resultState ? t.$LOCAL || "RESULT" == t.resultState ? t._e() : i("div", {
                    attrs: {
                        slot: "desc"
                    },
                    slot: "desc"
                }, [ i("v-container", {
                    attrs: {
                        "fill-height": "",
                        "grid-list-md": ""
                    }
                }, [ i("v-layout", {
                    attrs: {
                        "justify-center": ""
                    }
                }, [ i("v-flex", {
                    attrs: {
                        xs4: ""
                    }
                }, [ i("img", {
                    staticStyle: {
                        "max-width": "100%"
                    },
                    attrs: {
                        src: e(4833),
                        srcset: e(4833) + " 2x, " + e(5904) + " 3x"
                    }
                }) ]), t._v(" "), i("v-flex", {
                    attrs: {
                        xs8: ""
                    }
                }, [ i("h3", [ t._v("Error! ") ]), t._v(" "), i("p", [ t._v("Start a "), i("v-btn", {
                    attrs: {
                        to: "/search"
                    }
                }, [ t._v("New Search") ]), t._v("?") ], 1) ]) ], 1) ], 1) ], 1) : i("div", {
                    attrs: {
                        slot: "desc"
                    },
                    slot: "desc"
                }, [ i("v-container", {
                    attrs: {
                        "fill-height": "",
                        "grid-list-md": ""
                    }
                }, [ i("v-layout", {
                    attrs: {
                        "justify-center": ""
                    }
                }, [ i("v-flex", {
                    attrs: {
                        xs4: ""
                    }
                }, [ i("img", {
                    staticStyle: {
                        "max-width": "100%"
                    },
                    attrs: {
                        src: e(7969),
                        srcset: e(7969) + " 2x, " + e(5515) + " 3x"
                    }
                }) ]), t._v(" "), i("v-flex", {
                    attrs: {
                        xs8: ""
                    }
                }, [ i("h3", [ t._v("No hits found!") ]), t._v(" "), i("p", [ t._v("Start a "), i("v-btn", {
                    attrs: {
                        to: "/search"
                    }
                }, [ t._v("New Search") ]), t._v("?") ], 1) ]) ], 1) ], 1) ], 1) : i("div", {
                    attrs: {
                        slot: "desc"
                    },
                    slot: "desc"
                }, [ i("v-container", {
                    attrs: {
                        "fill-height": "",
                        "grid-list-md": ""
                    }
                }, [ i("v-layout", {
                    attrs: {
                        "justify-center": ""
                    }
                }, [ i("v-flex", {
                    attrs: {
                        xs4: ""
                    }
                }, [ i("img", {
                    staticStyle: {
                        "max-width": "100%"
                    },
                    attrs: {
                        src: e(4484),
                        srcset: e(4484) + " 2x, " + e(7940) + " 3x"
                    }
                }) ]), t._v(" "), i("v-flex", {
                    attrs: {
                        xs8: ""
                    }
                }, [ i("h3", [ t._v("Still Pending") ]), t._v(" "), i("p", [ t._v("Please wait a moment") ]) ]) ], 1) ], 1) ], 1), t._v(" "), t.hits && t.hits.results ? i("template", {
                    slot: "content"
                }, [ t.hits.results.length > 1 ? i("v-tabs", {
                    staticStyle: {
                        "margin-bottom": "2em"
                    },
                    attrs: {
                        color: t.selectedDatabases > 0 ? t.hits.results[t.selectedDatabases - 1].color : null,
                        "center-active": "",
                        grow: "",
                        "show-arrows": ""
                    },
                    on: {
                        change: function(n) {
                            return t.closeAlignment();
                        }
                    },
                    model: {
                        value: t.selectedDatabases,
                        callback: function(n) {
                            t.selectedDatabases = n;
                        },
                        expression: "selectedDatabases"
                    }
                }, [ i("v-tab", [ t._v("All databases") ]), t._v(" "), t._l(t.hits.results, (function(n) {
                    return i("v-tab", {
                        key: n.db
                    }, [ t._v(t._s(n.db) + " (" + t._s(n.alignments ? n.alignments.length : 0) + ")") ]);
                })) ], 2) : t._e(), t._v(" "), t._l(t.hits.results, (function(n, e) {
                    return 0 == t.selectedDatabases || e + 1 == t.selectedDatabases ? i("div", {
                        key: n.db
                    }, [ i("v-flex", {
                        staticClass: "d-flex",
                        style: {
                            "flex-direction": t.$vuetify.breakpoint.xsOnly ? "column" : null
                        }
                    }, [ i("h2", {
                        staticStyle: {
                            "margin-top": "0.5em",
                            "margin-bottom": "1em",
                            display: "inline-block"
                        }
                    }, [ i("span", {
                        staticStyle: {
                            "text-transform": "uppercase"
                        }
                    }, [ t._v(t._s(n.db)) ]), t._v(" "), i("small", [ t._v(t._s(n.alignments ? n.alignments.length : 0) + " hits") ]) ]), t._v(" "), i("v-btn-toggle", {
                        staticClass: "ml-auto",
                        attrs: {
                            mandatory: ""
                        },
                        model: {
                            value: t.tableMode,
                            callback: function(n) {
                                t.tableMode = n;
                            },
                            expression: "tableMode"
                        }
                    }, [ i("v-btn", [ t._v("\n                            Graphical\n                        ") ]), t._v(" "), i("v-btn", [ t._v("\n                            Numeric\n                        ") ]) ], 1) ], 1), t._v(" "), i("table", {
                        staticClass: "v-table result-table",
                        staticStyle: {
                            position: "relativ",
                            "margin-bottom": "3em"
                        }
                    }, [ i("thead", [ i("tr", [ i("th", {
                        class: "wide-" + (3 - n.hasDescription - n.hasTaxonomy)
                    }, [ t._v("Target") ]), t._v(" "), n.hasDescription ? i("th", {
                        staticClass: "wide-1"
                    }, [ t._v("\n                                Description\n                                "), i("v-tooltip", {
                        attrs: {
                            "open-delay": "300",
                            top: ""
                        },
                        scopedSlots: t._u([ {
                            key: "activator",
                            fn: function(n) {
                                var e = n.on;
                                return [ i("v-icon", t._g({
                                    staticStyle: {
                                        "font-size": "16px",
                                        float: "right"
                                    }
                                }, e), [ t._v(t._s(t.$MDI.HelpCircleOutline)) ]) ];
                            }
                        } ], null, !0)
                    }, [ t._v(" "), i("span", [ t._v("Triple click to select whole cell (for very long identifiers)") ]) ]) ], 1) : t._e(), t._v(" "), n.hasTaxonomy ? i("th", {
                        staticClass: "wide-1"
                    }, [ t._v("Scientific Name") ]) : t._e(), t._v(" "), i("th", {
                        staticClass: "thin"
                    }, [ t._v("Prob.") ]), t._v(" "), i("th", {
                        staticClass: "thin"
                    }, [ t._v("Seq. Id.") ]), t._v(" "), i("th", {
                        staticClass: "thin"
                    }, [ t._v(t._s("foldseek" == t.$APP && "tmalign" == t.mode ? "TM-score" : "E-Value")) ]), t._v(" "), 1 == t.tableMode ? i("th", {
                        staticClass: "thin"
                    }, [ t._v("Score") ]) : t._e(), t._v(" "), 1 == t.tableMode ? i("th", [ t._v("Query Pos.") ]) : t._e(), t._v(" "), 1 == t.tableMode ? i("th", [ t._v("Target Pos.") ]) : t._e(), t._v(" "), 0 == t.tableMode ? i("th", [ t._v("\n                                Position in query\n                                "), i("v-tooltip", {
                        attrs: {
                            "open-delay": "300",
                            top: ""
                        },
                        scopedSlots: t._u([ {
                            key: "activator",
                            fn: function(n) {
                                var e = n.on;
                                return [ i("v-icon", t._g({
                                    staticStyle: {
                                        "font-size": "16px",
                                        float: "right"
                                    }
                                }, e), [ t._v(t._s(t.$MDI.HelpCircleOutline)) ]) ];
                            }
                        } ], null, !0)
                    }, [ t._v(" "), i("span", [ t._v("The position of the aligned region of the target sequence in the query") ]) ]) ], 1) : t._e(), t._v(" "), i("th", {
                        staticClass: "alignment-action thin"
                    }, [ t._v("Alignment") ]) ]) ]), t._v(" "), i("tbody", t._l(n.alignments, (function(e, a) {
                        return i("tr", {
                            key: e.target + a,
                            class: [ "hit", {
                                active: e.active
                            } ]
                        }, [ i("td", {
                            staticClass: "long db",
                            style: "border-color: " + n.color,
                            attrs: {
                                "data-label": "Target"
                            }
                        }, [ i("a", {
                            staticClass: "anchor",
                            staticStyle: {
                                position: "absolute",
                                top: "0"
                            },
                            attrs: {
                                id: e.id
                            }
                        }), t._v(" "), i("a", {
                            attrs: {
                                href: e.href,
                                target: "_blank",
                                rel: "noopener",
                                title: e.target
                            }
                        }, [ t._v(t._s(e.target)) ]) ]), t._v(" "), n.hasDescription ? i("td", {
                            staticClass: "long",
                            attrs: {
                                "data-label": "Description"
                            }
                        }, [ i("span", {
                            attrs: {
                                title: e.description
                            }
                        }, [ t._v(t._s(e.description)) ]) ]) : t._e(), t._v(" "), n.hasTaxonomy ? i("td", {
                            staticClass: "long",
                            attrs: {
                                "data-label": "Taxonomy"
                            }
                        }, [ i("a", {
                            attrs: {
                                href: "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=" + e.taxId,
                                target: "_blank",
                                rel: "noopener",
                                title: e.taxName
                            }
                        }, [ t._v(t._s(e.taxName)) ]) ]) : t._e(), t._v(" "), i("td", {
                            staticClass: "thin",
                            attrs: {
                                "data-label": "Probability"
                            }
                        }, [ t._v(t._s(e.prob)) ]), t._v(" "), i("td", {
                            staticClass: "thin",
                            attrs: {
                                "data-label": "Sequence Identity"
                            }
                        }, [ t._v(t._s(e.seqId)) ]), t._v(" "), i("td", {
                            staticClass: "thin",
                            attrs: {
                                "data-label": "foldseek" == t.$APP && "tmalign" == t.mode ? "TM-score" : "E-Value"
                            }
                        }, [ t._v(t._s(e.eval)) ]), t._v(" "), 1 == t.tableMode ? i("td", {
                            staticClass: "thin",
                            attrs: {
                                "data-label": "Score"
                            }
                        }, [ t._v(t._s(e.score)) ]) : t._e(), t._v(" "), 1 == t.tableMode ? i("td", {
                            staticClass: "thin",
                            attrs: {
                                "data-label": "Query Position"
                            }
                        }, [ t._v(t._s(e.qStartPos) + "-" + t._s(e.qEndPos) + " (" + t._s(e.qLen) + ")") ]) : t._e(), t._v(" "), 1 == t.tableMode ? i("td", {
                            staticClass: "thin",
                            attrs: {
                                "data-label": "Target Position"
                            }
                        }, [ t._v(t._s(e.dbStartPos) + "-" + t._s(e.dbEndPos) + " (" + t._s(e.dbLen) + ")") ]) : t._e(), t._v(" "), 0 == t.tableMode ? i("td", {
                            staticClass: "graphical",
                            attrs: {
                                "data-label": "Position"
                            }
                        }, [ i("Ruler", {
                            attrs: {
                                length: e.qLen,
                                start: e.qStartPos,
                                end: e.qEndPos,
                                color: e.color,
                                label: 0 == a
                            }
                        }) ], 1) : t._e(), t._v(" "), i("td", {
                            staticClass: "alignment-action thin"
                        }, [ i("button", {
                            staticClass: "v-btn v-btn--icon v-btn--round v-btn--text v-size--default",
                            class: {
                                "v-btn--outlined": t.alignment && e.target == t.alignment.target,
                                "theme--dark": t.$vuetify.theme.dark
                            },
                            attrs: {
                                type: "button"
                            },
                            on: {
                                click: function(n) {
                                    return t.showAlignment(e, n);
                                }
                            }
                        }, [ i("span", {
                            staticClass: "v-btn__content"
                        }, [ i("span", {
                            staticClass: "v-icon notranslate theme--dark",
                            attrs: {
                                "aria-hidden": "true"
                            }
                        }, [ i("svg", {
                            staticClass: "v-icon__svg",
                            attrs: {
                                xmlns: "http://www.w3.org/2000/svg",
                                viewBox: "0 0 24 24",
                                role: "img",
                                "aria-hidden": "true"
                            }
                        }, [ i("path", {
                            attrs: {
                                d: "M5,13H19V11H5M3,17H17V15H3M7,7V9H21V7"
                            }
                        }) ]) ]) ]) ]) ]) ]);
                    })), 0) ]) ], 1) : t._e();
                })) ], 2) : t._e() ], 2) ], 1) ], 1), t._v(" "), i("portal", [ null != t.alignment ? i("panel", {
                    staticClass: "alignment",
                    style: "top: " + t.alnBoxOffset + "px"
                }, [ i("AlignmentPanel", {
                    key: "ap-" + t.alignment.id,
                    attrs: {
                        slot: "content",
                        alignment: t.alignment,
                        lineLen: t.fluidLineLen,
                        hits: t.hits
                    },
                    slot: "content"
                }) ], 1) : t._e() ], 1) ], 1);
            };
            x._withStripped = !0;
            var I = e(917), S = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("div", {
                    staticClass: "alignment-wrapper-outer",
                    attrs: {
                        slot: "content"
                    },
                    slot: "content"
                }, [ e("Alignment", {
                    key: "aln2-" + t.alignment.id,
                    attrs: {
                        alignment: t.alignment,
                        lineLen: t.lineLen,
                        queryMap: t.queryMap,
                        targetMap: t.targetMap
                    },
                    on: {
                        selected: t.setUserSelection
                    }
                }), t._v(" "), "foldseek" == t.$APP ? e("div", {
                    staticClass: "alignment-structure-wrapper"
                }, [ e("StructureViewer", {
                    key: "struc2-" + t.alignment.id,
                    ref: "structureViewer",
                    attrs: {
                        alignment: t.alignment,
                        queryMap: t.queryMap,
                        targetMap: t.targetMap,
                        hits: t.hits,
                        bgColorLight: "white",
                        bgColorDark: "#1E1E1E",
                        qColor: "lightgrey",
                        tColor: "red",
                        qRepr: "cartoon",
                        tRepr: "cartoon"
                    }
                }) ], 1) : t._e() ], 1);
            };
            S._withStripped = !0;
            var T = e(8152), N = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("div", {
                    staticClass: "alignment-wrapper-inner"
                }, [ t._l(Math.max(1, Math.ceil(t.alignment.alnLength / t.lineLen)), (function(n) {
                    return e("span", {
                        key: n,
                        staticClass: "monospace"
                    }, [ e("span", {
                        staticClass: "line"
                    }, [ t._v("\n            Q " + t._s(t.padNumber(t.getQueryRowStartPos(n), (Math.max(t.alignment.qStartPos, t.alignment.dbStartPos) + t.alignment.alnLength + "").length, " ")) + " "), e("span", {
                        staticClass: "residues"
                    }, [ t._v(t._s(t.alignment.qAln.substring((n - 1) * t.lineLen, (n - 1) * t.lineLen + t.lineLen))) ]), t._v(" "), e("br"), t._v("\n            " + t._s(" ".repeat(3 + (Math.max(t.alignment.qStartPos, t.alignment.dbStartPos) + t.alignment.alnLength + "").length))), e("span", {
                        staticClass: "residues"
                    }, [ t._v(t._s(t.formatAlnDiff(t.alignment.qAln.substring((n - 1) * t.lineLen, (n - 1) * t.lineLen + t.lineLen), t.alignment.dbAln.substring((n - 1) * t.lineLen, (n - 1) * t.lineLen + t.lineLen)))) ]), t._v(" "), e("br"), t._v("\n            T " + t._s(t.padNumber(t.getTargetRowStartPos(n), (Math.max(t.alignment.qStartPos, t.alignment.dbStartPos) + t.alignment.alnLength + "").length, " ")) + " "), e("span", {
                        staticClass: "residues",
                        on: {
                            pointerup: function(e) {
                                return t.onSelectText(n);
                            }
                        }
                    }, [ t._v(t._s(t.alignment.dbAln.substring((n - 1) * t.lineLen, (n - 1) * t.lineLen + t.lineLen))) ]) ]), e("br") ]);
                })), t._v(" "), "foldseek" == t.$APP ? e("small", {
                    staticStyle: {
                        float: "right"
                    }
                }, [ t._v("Select target residues to highlight their structure") ]) : t._e() ], 2);
            };
            N._withStripped = !0;
            var L = [ "AG", "AS", "DE", "DN", "ED", "EK", "EQ", "FL", "FM", "FW", "FY", "GA", "HN", "HQ", "HY", "IL", "IM", "IV", "KE", "KQ", "KR", "LF", "LI", "LM", "LV", "MF", "MI", "ML", "MV", "ND", "NH", "NQ", "NS", "QE", "QH", "QK", "QN", "QR", "RK", "RQ", "SA", "SN", "ST", "TS", "VI", "VL", "VM", "WF", "WY", "YF", "YH", "YW" ];
            const E = {
                props: [ "alignment", "lineLen", "queryMap", "targetMap" ],
                methods: {
                    getQueryIndex: function(t) {
                        return this.queryMap[t];
                    },
                    getTargetIndex: function(t) {
                        return this.targetMap[t];
                    },
                    getFirstResidueNumber: function(t, n) {
                        for (var e = this.lineLen * (n - 1); null === t[e]; ) e--;
                        return t[e];
                    },
                    getQueryRowStartPos: function(t) {
                        return this.getFirstResidueNumber(this.queryMap, t);
                    },
                    getTargetRowStartPos: function(t) {
                        return this.getFirstResidueNumber(this.targetMap, t);
                    },
                    formatAlnDiff: function(t, n) {
                        if (t.length != n.length) return "";
                        for (var e = "", i = 0; i < t.length; i++) t[i] == n[i] ? e += t[i] : -1 != L.indexOf(t[i] + n[i]) ? e += "+" : e += " ";
                        return e;
                    },
                    padNumber: function(t, n, e) {
                        return Array(n - String(t).length + 1).join(e || "0") + t;
                    },
                    onSelectText: function(t) {
                        var n = window.getSelection(), e = [ n.anchorOffset, n.focusOffset ].sort((function(t, n) {
                            return t - n;
                        })), i = (0, T.Z)(e, 2), a = i[0], r = i[1] - a, s = (t - 1) * this.lineLen + a, o = s + r - 1, l = function(t, n, e) {
                            for (var i = null, a = null, r = n; r <= e; r++) {
                                var s = t[r];
                                null !== s && (null === i && (i = s), a = s);
                            }
                            return [ i, a ];
                        }(this.targetMap, s, o), A = (0, T.Z)(l, 2), c = A[0], d = A[1];
                        this.$emit("selected", [ c, d ]);
                    }
                }
            };
            e(603);
            var D = (0, M.Z)(E, N, [], !1, null, null, null);
            D.options.__file = "frontend/Alignment.vue";
            function O(t, n) {
                for (var e = Array(n.length), i = 0, a = 0; i < n.length; i++) "-" === n[i] ? (e[i] = null, 
                a++) : e[i] = t + i - a;
                return e;
            }
            const R = {
                components: {
                    StructureViewer: function() {
                        return Promise.resolve().then(e.bind(e, 8992));
                    },
                    Alignment: D.exports
                },
                data: function() {
                    return {
                        queryMap: null,
                        targetMap: null
                    };
                },
                props: {
                    alignment: {
                        type: Object,
                        required: !0
                    },
                    lineLen: {
                        type: Number,
                        required: !0
                    },
                    hits: {
                        type: Object
                    }
                },
                methods: {
                    setUserSelection: function(t) {
                        var n = (0, T.Z)(t, 2), e = n[0], i = n[1];
                        this.alignment && this.$refs.structureViewer.setSelectionData(e, i);
                    },
                    updateMaps: function() {
                        this.alignment && (this.queryMap = O(this.alignment.qStartPos, this.alignment.qAln), 
                        this.targetMap = O(this.alignment.dbStartPos, this.alignment.dbAln));
                    }
                },
                watch: {
                    alignment: function() {
                        this.updateMaps();
                    }
                },
                beforeMount: function() {
                    this.updateMaps();
                }
            };
            e(2530);
            var k = (0, M.Z)(R, S, [], !1, null, null, null);
            k.options.__file = "frontend/AlignmentPanel.vue";
            const B = k.exports;
            var P = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("div", {
                    staticClass: "ruler"
                }, [ e("div", {
                    staticClass: "query",
                    class: {
                        reversed: t.reversed
                    },
                    style: {
                        left: t.queryLeft + "%",
                        right: t.queryRight + "%"
                    }
                }, [ e("div", {
                    staticClass: "chevron-start",
                    style: {
                        "background-color": t.color
                    }
                }), t._v(" "), e("div", {
                    staticClass: "chevron-mid",
                    style: {
                        "background-color": t.color
                    }
                }), t._v(" "), e("div", {
                    staticClass: "chevron-end",
                    style: {
                        "background-color": t.color
                    }
                }) ]), t._v(" "), e("div", {
                    staticClass: "tick-label",
                    style: {
                        left: t.queryLeft + "%"
                    }
                }, [ t._v(t._s(t.minStart)) ]), t._v(" "), e("div", {
                    staticClass: "tick-label",
                    style: {
                        right: t.queryRight + "%",
                        "margin-left": 0,
                        "margin-right": "-25px"
                    }
                }, [ t._v(t._s(t.maxEnd)) ]) ]);
            };
            P._withStripped = !0;
            const z = {
                props: {
                    length: Number,
                    start: Number,
                    end: Number,
                    color: String,
                    label: Boolean,
                    tickInterval: {
                        type: Number,
                        default: 10
                    }
                },
                computed: {
                    minStart: function() {
                        return Math.min(this.start, this.end);
                    },
                    maxEnd: function() {
                        return Math.max(this.start, this.end);
                    },
                    reversed: function() {
                        return this.start > this.end;
                    },
                    queryLeft: function() {
                        return (this.minStart - 1) / this.length * 100;
                    },
                    queryRight: function() {
                        return 100 - this.maxEnd / this.length * 100;
                    },
                    numTicks: function() {
                        return 3;
                    },
                    ticks: function() {
                        var t = this;
                        return Array.from({
                            length: this.numTicks + 1
                        }, (function(n, e) {
                            return e / t.numTicks * 100;
                        }));
                    }
                }
            };
            e(5941);
            var Z = (0, M.Z)(z, P, [], !1, null, "2b7861b2", null);
            Z.options.__file = "frontend/Ruler.vue";
            const _ = Z.exports;
            function U(t, n, e) {
                var i;
                return function() {
                    var a = this, r = arguments, s = e && !i;
                    clearTimeout(i), i = setTimeout((function() {
                        i = null, e || t.apply(a, r);
                    }), n), s && t.apply(a, r);
                };
            }
            function j(t) {
                for (var n = 0; t; ) n += t.offsetTop, t = t.offsetParent;
                return n;
            }
            const q = {
                name: "result",
                components: {
                    Panel: I.Z,
                    AlignmentPanel: B,
                    Ruler: _
                },
                data: function() {
                    return {
                        alignment: null,
                        activeTarget: null,
                        alnBoxOffset: 0,
                        selectedDatabases: 0,
                        tableMode: 0
                    };
                },
                props: {
                    ticket: "",
                    error: "",
                    mode: "",
                    hits: null
                },
                created: function() {
                    window.addEventListener("resize", this.handleAlignmentBoxResize, {
                        passive: !0
                    });
                },
                beforeDestroy: function() {
                    window.removeEventListener("resize", this.handleAlignmentBoxResize);
                },
                computed: {
                    fluidLineLen: function() {
                        return this.$vuetify.breakpoint.xsOnly ? 30 : this.$vuetify.breakpoint.smAndDown ? 40 : 80;
                    },
                    filteredResults: function() {
                        return this.hits ? 0 === this.selectedDatabases ? this.hits.results : [ this.hits.results[this.selectedDatabases - 1] ] : [];
                    },
                    resultState: function() {
                        if (null == this.hits && "" == this.error) return "PENDING";
                        if (!this.hits.results) return "ERROR";
                        if (0 == this.hits.results.length) return "EMPTY";
                        for (var t in this.hits.results) if (null != this.hits.results[t].alignments) return "RESULT";
                        return "ERROR";
                    }
                },
                methods: {
                    showAlignment: function(t, n) {
                        this.alignment === t ? this.closeAlignment() : (this.alignment = t, this.activeTarget = n.target.closest(".hit"), 
                        this.alnBoxOffset = j(this.activeTarget) + this.activeTarget.offsetHeight);
                    },
                    closeAlignment: function() {
                        this.alignment = null, this.activeTarget = null;
                    },
                    handleAlignmentBoxResize: U((function() {
                        null != this.activeTarget && (this.alnBoxOffset = j(this.activeTarget) + this.activeTarget.offsetHeight);
                    }), 32, !1)
                }
            };
            e(5264);
            var G = e(3453), Q = e.n(G), V = e(5934), F = e(6584), H = e(6530), Y = e(683), W = e(4786), J = e(9456), K = e(756), X = e(7849), $ = e(1562), tt = (0, 
            M.Z)(q, x, [], !1, null, null, null);
            Q()(tt, {
                VBtn: V.Z,
                VBtnToggle: F.Z,
                VContainer: H.Z,
                VFlex: Y.Z,
                VIcon: W.Z,
                VLayout: J.Z,
                VTab: K.Z,
                VTabs: X.Z,
                VTooltip: $.Z
            }), tt.options.__file = "frontend/ResultView.vue";
            const nt = tt.exports;
            var et = function() {
                var t = this, n = t.$createElement, i = t._self._c || n;
                return i("div", [ t.$LOCAL ? t._e() : i("v-navigation-drawer", {
                    ref: "drawer",
                    attrs: {
                        stateless: "",
                        app: "",
                        permanent: "",
                        clipped: "",
                        "mini-variant": t.mini,
                        "expand-on-hover": !1
                    }
                }, [ t.$LOCAL ? t._e() : i("v-list", [ i("v-list-item", {
                    attrs: {
                        to: "/search"
                    }
                }, [ i("v-list-item-action", [ i("v-icon", [ t._v(t._s(t.$MDI.Magnify)) ]) ], 1), t._v(" "), i("v-list-item-content", [ i("v-list-item-title", [ t._v("Search") ]) ], 1) ], 1), t._v(" "), i("v-list-item", [ i("v-list-item-action", [ i("v-icon", [ t._v(t._s(t.$MDI.Upload)) ]) ], 1), t._v(" "), i("v-list-item-content", [ i("v-list-item-title", [ t._v("Upload JSON") ]) ], 1) ], 1), t._v(" "), i("router-view", {
                    attrs: {
                        name: "sidebar"
                    }
                }), t._v(" "), t.$LOCAL ? t._e() : i("history"), t._v(" "), t.$ELECTRON ? i("v-list-item", {
                    attrs: {
                        to: "/preferences"
                    }
                }, [ i("v-list-item-action", [ i("v-icon", [ t._v(t._s(t.$MDI.Tune)) ]) ], 1), t._v(" "), i("v-list-item-content", [ i("v-list-item-title", [ t._v("Preferences") ]) ], 1) ], 1) : t._e() ], 1) ], 1), t._v(" "), i("v-app-bar", {
                    class: [ "ml-0", "pl-3", t.$ELECTRON ? "pt-2" : null ],
                    style: {
                        "-webkit-app-region": t.$ELECTRON ? "drag" : null,
                        "-webkit-user-select": t.$ELECTRON ? "none" : null
                    },
                    attrs: {
                        app: "",
                        height: t.$ELECTRON ? "72px" : "48px",
                        fixed: "",
                        "clipped-left": ""
                    },
                    nativeOn: {
                        dblclick: function(n) {
                            return t.electronHandleTitleBarDoubleClick();
                        }
                    }
                }, [ t.$LOCAL ? t._e() : i("v-app-bar-nav-icon", {
                    attrs: {
                        "input-value": t.mini ? void 0 : "activated"
                    },
                    on: {
                        click: function(n) {
                            return n.stopPropagation(), t.toggleMini.apply(null, arguments);
                        }
                    }
                }), t._v(" "), i("v-app-bar-title", [ t.$LOCAL ? t._e() : i("router-link", {
                    staticStyle: {
                        color: "inherit",
                        "text-decoration": "none"
                    },
                    attrs: {
                        to: "/"
                    }
                }, [ t._v(t._s(t.$STRINGS.APP_NAME) + " Search") ]), t._v(" "), t.$LOCAL ? i("span", [ t._v(t._s(t.$STRINGS.APP_NAME) + " Search") ]) : t._e() ], 1), t._v(" "), "mmseqs" == t.$APP ? i("object", {
                    staticStyle: {
                        "margin-left": "8px",
                        display: "inline-block",
                        width: "38px",
                        height: "38px",
                        "vertical-align": "middle"
                    },
                    attrs: {
                        type: "image/svg+xml",
                        data: e(2881),
                        "aria-hidden": "true"
                    }
                }, [ i("img", {
                    staticStyle: {
                        "max-width": "100%"
                    },
                    attrs: {
                        src: e(7018)
                    }
                }) ]) : t._e(), t._v(" "), "foldseek" == t.$APP ? i("img", {
                    staticStyle: {
                        "margin-left": "8px",
                        display: "inline-block",
                        width: "48px",
                        height: "48px",
                        "vertical-align": "middle"
                    },
                    attrs: {
                        src: e(6617),
                        "aria-hidden": "true"
                    }
                }) : t._e(), t._v(" "), i("v-spacer"), t._v(" "), t.$ELECTRON ? t._e() : t._m(0) ], 1) ], 1);
            };
            et._withStripped = !0;
            var it = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return t.items && t.items.length > 0 ? e("v-list-group", {
                    attrs: {
                        "no-action": "",
                        "prepend-icon": t.$MDI.History
                    },
                    model: {
                        value: t.drawer,
                        callback: function(n) {
                            t.drawer = n;
                        },
                        expression: "drawer"
                    }
                }, [ e("template", {
                    slot: "activator"
                }, [ e("v-list-item-content", [ e("v-list-item-title", [ t._v("\n                History\n            ") ]), t._v(" "), t.drawer ? e("v-list-item-subtitle", {
                    staticClass: "ml-n1",
                    on: {
                        click: function(t) {
                            t.preventDefault();
                        }
                    }
                }, [ e("button", {
                    style: {
                        opacity: 0 == t.page ? .6 : 1
                    },
                    on: {
                        click: function(n) {
                            return n.preventDefault(), t.previous();
                        }
                    }
                }, [ e("v-icon", {
                    staticStyle: {
                        transform: "inherit"
                    }
                }, [ t._v(t._s(t.$MDI.ChevronLeft)) ]) ], 1), t._v(" "), e("button", {
                    style: {
                        opacity: (t.page + 1) * t.limit >= t.items.length ? .6 : 1
                    },
                    on: {
                        click: function(n) {
                            return n.preventDefault(), t.next();
                        }
                    }
                }, [ e("v-icon", {
                    staticStyle: {
                        transform: "inherit"
                    }
                }, [ t._v(t._s(t.$MDI.ChevronRight)) ]) ], 1) ]) : t._e() ], 1) ], 1), t._v(" "), t._l(t.items.slice(t.page * t.limit, (t.page + 1) * t.limit), (function(n, i) {
                    return e("v-list-item", {
                        key: i,
                        class: {
                            "list__item--highlighted": n.id == t.current
                        },
                        staticStyle: {
                            "padding-left": "16px"
                        },
                        attrs: {
                            to: t.formattedRoute(n)
                        }
                    }, [ e("v-list-item-icon", [ "COMPLETE" == n.status ? e("identicon", {
                        attrs: {
                            hash: n.id
                        }
                    }) : "RUNNING" == n.status || "PENDING" == n.status ? e("v-icon", {
                        attrs: {
                            large: ""
                        }
                    }, [ t._v(t._s(t.$MDI.ClockOutline)) ]) : (n.status, e("v-icon", {
                        attrs: {
                            large: ""
                        }
                    }, [ t._v(t._s(t.$MDI.HelpCircleOutline)) ])) ], 1), t._v(" "), e("v-list-item-content", [ e("v-list-item-title", [ t._v("\n                " + t._s(t.formattedDate(n.time)) + "\n            ") ]), t._v(" "), e("v-list-item-subtitle", [ e("span", {
                        staticClass: "mono"
                    }, [ t._v(t._s(n.id)) ]) ]) ], 1) ], 1);
                })) ], 2) : t._e();
            };
            it._withStripped = !0;
            var at = function() {
                var t = this, n = t.$createElement;
                return (t._self._c || n)("img", {
                    attrs: {
                        src: "data:image/svg+xml;base64," + t.makeData(t.hash, t.size),
                        width: t.size,
                        height: t.size
                    }
                });
            };
            at._withStripped = !0;
            const rt = function() {
                var t = function(t, n) {
                    if ("string" != typeof t || t.length < 15) throw "A hash of at least 15 characters is required.";
                    this.defaults = {
                        background: [ 240, 240, 240, 255 ],
                        margin: .08,
                        size: 64,
                        saturation: .7,
                        brightness: .5,
                        format: "png"
                    }, this.options = "object" === (0, r.Z)(n) ? n : this.defaults, "number" == typeof arguments[1] && (this.options.size = arguments[1]), 
                    arguments[2] && (this.options.margin = arguments[2]), this.hash = t, this.background = this.options.background || this.defaults.background, 
                    this.size = this.options.size || this.defaults.size, this.format = this.options.format || this.defaults.format, 
                    this.margin = void 0 !== this.options.margin ? this.options.margin : this.defaults.margin;
                    var e = parseInt(this.hash.substr(-7), 16) / 268435455, i = this.options.saturation || this.defaults.saturation, a = this.options.brightness || this.defaults.brightness;
                    this.foreground = this.options.foreground || this.hsl2rgb(e, i, a);
                };
                t.prototype = {
                    background: null,
                    foreground: null,
                    hash: null,
                    margin: null,
                    size: null,
                    format: null,
                    image: function() {
                        return this.isSvg() ? new n(this.size, this.foreground, this.background) : new null(this.size, this.size, 256);
                    },
                    render: function() {
                        var t, n, e = this.image(), i = this.size, a = Math.floor(i * this.margin), r = Math.floor((i - 2 * a) / 5), s = Math.floor((i - 5 * r) / 2), o = e.color.apply(e, this.background), l = e.color.apply(e, this.foreground);
                        for (t = 0; t < 15; t++) n = parseInt(this.hash.charAt(t), 16) % 2 ? o : l, t < 5 ? this.rectangle(2 * r + s, t * r + s, r, r, n, e) : t < 10 ? (this.rectangle(1 * r + s, (t - 5) * r + s, r, r, n, e), 
                        this.rectangle(3 * r + s, (t - 5) * r + s, r, r, n, e)) : t < 15 && (this.rectangle(0 * r + s, (t - 10) * r + s, r, r, n, e), 
                        this.rectangle(4 * r + s, (t - 10) * r + s, r, r, n, e));
                        return e;
                    },
                    rectangle: function(t, n, e, i, a, r) {
                        var s, o;
                        if (this.isSvg()) r.rectangles.push({
                            x: t,
                            y: n,
                            w: e,
                            h: i,
                            color: a
                        }); else for (s = t; s < t + e; s++) for (o = n; o < n + i; o++) r.buffer[r.index(s, o)] = a;
                    },
                    hsl2rgb: function(t, n, e) {
                        return [ 255 * (n = [ e += n *= e < .5 ? e : 1 - e, e - (t *= 6) % 1 * n * 2, e -= n *= 2, e, e + t % 1 * n, e + n ])[~~t % 6], 255 * n[(16 | t) % 6], 255 * n[(8 | t) % 6] ];
                    },
                    toString: function(t) {
                        return t ? this.render().getDump() : this.render().getBase64();
                    },
                    isSvg: function() {
                        return this.format.match(/svg/i);
                    }
                };
                var n = function(t, n, e) {
                    this.size = t, this.foreground = this.color.apply(this, n), this.background = this.color.apply(this, e), 
                    this.rectangles = [];
                };
                return n.prototype = {
                    size: null,
                    foreground: null,
                    background: null,
                    rectangles: null,
                    color: function(t, n, e, i) {
                        var a = [ t, n, e ].map(Math.round);
                        return a.push(i >= 0 && i <= 255 ? i / 255 : 1), "rgba(" + a.join(",") + ")";
                    },
                    getDump: function() {
                        var t, n, e, i = this.foreground, a = this.background, r = .005 * this.size;
                        for (n = "<svg xmlns='http://www.w3.org/2000/svg' width='" + this.size + "' height='" + this.size + "' style='background-color:" + a + ";'><g style='fill:" + i + "; stroke:" + i + "; stroke-width:" + r + ";'>", 
                        t = 0; t < this.rectangles.length; t++) (e = this.rectangles[t]).color != a && (n += "<rect  x='" + e.x + "' y='" + e.y + "' width='" + e.w + "' height='" + e.h + "'/>");
                        return n += "</g></svg>";
                    },
                    getBase64: function() {
                        return btoa(this.getDump());
                    }
                }, t;
            }();
            const st = {
                name: "identicon",
                props: {
                    hash: {
                        default: "",
                        type: String
                    },
                    size: {
                        default: 32,
                        type: Number
                    }
                },
                methods: {
                    makeData: function(t, n) {
                        return new rt(function(t) {
                            for (var n = 0, e = 0; e < t.length; ++e) n = 31 * n + t[e].charCodeAt(0);
                            return n.toString(16).slice(0, 14) + "" + n.toString(16)[0];
                        }(t), {
                            background: [ 0, 0, 0, 0 ],
                            margin: 0,
                            size: n,
                            format: "svg"
                        }).toString();
                    }
                }
            }, ot = st;
            var lt = (0, M.Z)(ot, at, [], !1, null, null, null);
            lt.options.__file = "frontend/Identicon.vue";
            const At = lt.exports;
            var ct = !1;
            try {
                void 0 !== window.localStorage && (ct = !0);
            } catch (t) {}
            const dt = {
                components: {
                    Identicon: At
                },
                data: function() {
                    return {
                        current: "",
                        drawer: !1,
                        error: !1,
                        items: [],
                        page: 0,
                        limit: 7
                    };
                },
                mounted: function() {},
                created: function() {
                    this.fetchData();
                },
                watch: {
                    $route: function(t, n) {
                        n.path != t.path && this.fetchData();
                    },
                    items: function(t) {
                        ct && (localStorage.history = JSON.stringify(t));
                    },
                    drawer: function(t, n) {
                        1 == t && this.$root.$emit("multi", !0);
                    }
                },
                methods: {
                    previous: function() {
                        0 != this.page && (this.page -= 1);
                    },
                    next: function() {
                        (this.page + 1) * this.limit > this.items.length || (this.page += 1);
                    },
                    fetchData: U((function() {
                        var t, n = this;
                        this.current = this.$route.params.ticket, this.error = !1, t = ct && localStorage.history ? JSON.parse(localStorage.history) : [];
                        var e = [], i = !1;
                        for (var a in t) this.current == t[a].id && (i = !0), e.push(t[a].id), t[a].status = "UNKNOWN";
                        null != this.current && 0 == i && (e.unshift(this.current), t.unshift({
                            id: this.current,
                            status: "UNKNOWN",
                            time: +new Date
                        })), this.$axios.post("api/tickets", function(t) {
                            for (var n = new URLSearchParams(t), e = function() {
                                var t = (0, T.Z)(a[i], 2), e = t[0], r = t[1];
                                Array.isArray(r) && (n.delete(e), r.forEach((function(t) {
                                    return n.append(e + "[]", t);
                                })));
                            }, i = 0, a = Object.entries(t); i < a.length; i++) e();
                            return n.toString();
                        }({
                            tickets: e
                        })).then((function(e) {
                            var i = e.data, a = +new Date, r = [], s = !1;
                            for (var o in i) {
                                var l = !1;
                                if ("COMPLETE" == i[o].status ? l = !0 : "UNKNOWN" == i[o].status ? l = !1 : a - t[o].time < 6048e5 && (l = !0), 
                                "PENDING" != i[o].status && "RUNNING" != i[o].status || (s = !0), l) {
                                    var A = t[o];
                                    A.status = i[o].status, r.push(A);
                                }
                            }
                            n.items = r, s && setTimeout(n.fetchData.bind(n), 5e3);
                        }), (function() {
                            n.error = !0;
                        }));
                    }), 16, !0),
                    formattedRoute: function(t) {
                        return "COMPLETE" == t.status ? "/result/" + t.id + "/0" : "/queue/" + t.id;
                    },
                    formattedDate: function(t) {
                        var n = new Date(t), e = n.getMonth() + 1, i = n.getDate(), a = n.getHours(), r = n.getMinutes();
                        return e = (e < 10 ? "0" : "") + e, i = (i < 10 ? "0" : "") + i, a = (a < 10 ? "0" : "") + a, 
                        r = (r < 10 ? "0" : "") + r, n.getFullYear() + "-" + e + "-" + i + " " + a + ":" + r;
                    }
                }
            };
            var ut = e(3308), ht = e(3347), pt = e(9623), gt = e(3560), mt = (0, M.Z)(dt, it, [], !1, null, null, null);
            Q()(mt, {
                VIcon: W.Z,
                VListGroup: ut.Z,
                VListItem: ht.Z,
                VListItemContent: pt.km,
                VListItemIcon: gt.Z,
                VListItemSubtitle: pt.oZ,
                VListItemTitle: pt.V9
            }), mt.options.__file = "frontend/History.vue";
            const vt = {
                components: {
                    History: mt.exports
                },
                data: function() {
                    return {
                        mini: !0
                    };
                },
                created: function() {
                    this.$root.$on("multi", this.shouldExpand);
                },
                mounted: function() {
                    console.log(this);
                },
                beforeDestroy: function() {
                    this.$root.$off("multi", this.shouldExpand);
                },
                methods: {
                    shouldExpand: function(t) {
                        this.mini = !t;
                    },
                    toggleMini: function() {
                        this.mini = !this.mini;
                    },
                    electronHandleTitleBarDoubleClick: function() {
                        this.handleTitleBarDoubleClick();
                    }
                }
            };
            e(4449);
            var ft = e(9085), bt = e(5078), Ct = e(8895), Mt = e(2545), yt = e(3444), wt = e(9681), xt = e(2515), It = e(3845), St = (0, 
            M.Z)(vt, et, [ function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("v-toolbar-items", {
                    staticClass: "hidden-sm-and-down"
                }, t._l(t.$STRINGS.NAV_URL_COUNT - 0, (function(n) {
                    return e("v-btn", {
                        key: n,
                        attrs: {
                            text: "",
                            rel: "external noopener",
                            target: "_blank",
                            href: t.$STRINGS["NAV_URL_" + n]
                        }
                    }, [ t._v(t._s(t.$STRINGS["NAV_TITLE_" + n])) ]);
                })), 1);
            } ], !1, null, "5976e89a", null);
            Q()(St, {
                VAppBar: ft.Z,
                VAppBarNavIcon: bt.Z,
                VAppBarTitle: Ct.Z,
                VBtn: V.Z,
                VIcon: W.Z,
                VList: Mt.Z,
                VListItem: ht.Z,
                VListItemAction: yt.Z,
                VListItemContent: pt.km,
                VListItemTitle: pt.V9,
                VNavigationDrawer: wt.Z,
                VSpacer: xt.Z,
                VToolbarItems: It.lj
            }), St.options.__file = "frontend/Navigation.vue";
            function Tt(t, n) {
                var e = "undefined" != typeof Symbol && t[Symbol.iterator] || t["@@iterator"];
                if (!e) {
                    if (Array.isArray(t) || (e = function(t, n) {
                        if (!t) return;
                        if ("string" == typeof t) return Nt(t, n);
                        var e = Object.prototype.toString.call(t).slice(8, -1);
                        "Object" === e && t.constructor && (e = t.constructor.name);
                        if ("Map" === e || "Set" === e) return Array.from(t);
                        if ("Arguments" === e || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(e)) return Nt(t, n);
                    }(t)) || n && t && "number" == typeof t.length) {
                        e && (t = e);
                        var i = 0, a = function() {};
                        return {
                            s: a,
                            n: function() {
                                return i >= t.length ? {
                                    done: !0
                                } : {
                                    done: !1,
                                    value: t[i++]
                                };
                            },
                            e: function(t) {
                                throw t;
                            },
                            f: a
                        };
                    }
                    throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.");
                }
                var r, s = !0, o = !1;
                return {
                    s: function() {
                        e = e.call(t);
                    },
                    n: function() {
                        var t = e.next();
                        return s = t.done, t;
                    },
                    e: function(t) {
                        o = !0, r = t;
                    },
                    f: function() {
                        try {
                            s || null == e.return || e.return();
                        } finally {
                            if (o) throw r;
                        }
                    }
                };
            }
            function Nt(t, n) {
                (null == n || n > t.length) && (n = t.length);
                for (var e = 0, i = new Array(n); e < n; e++) i[e] = t[e];
                return i;
            }
            const Lt = {
                name: "result",
                mixins: [ w ],
                components: {
                    ResultView: nt,
                    Navigation: St.exports
                },
                data: function() {
                    return {
                        currentIndex: 0
                    };
                },
                mounted: function() {
                    var t = this;
                    document.onreadystatechange = function() {
                        if ("complete" == document.readyState) {
                            var n = document.getElementById("data");
                            if (!n) return null;
                            var e = JSON.parse(n.textContent);
                            t.fetchData(e);
                        }
                    };
                },
                computed: {
                    currentResult: function() {
                        return null === this.hits ? null : this.hits[this.currentIndex];
                    },
                    currentQuery: function() {
                        return null === this.hits ? "" : this.hits[this.currentIndex].query.header;
                    }
                },
                methods: {
                    changeResult: function(t) {
                        this.currentIndex = t, this.setColorScheme();
                    },
                    uploadData: function(t) {
                        var n = this;
                        if (t) {
                            var e = new FileReader;
                            e.addEventListener("load", (function(t) {
                                var e = JSON.parse(t.target.result);
                                n.fetchData(e);
                            })), e.readAsText(t);
                        }
                    },
                    downloadData: function() {
                        if (!this.hits) return null;
                        var t = JSON.stringify(this.hits), n = new Blob([ t ], {
                            type: "application/json"
                        }), e = document.createElement("a"), i = (new Date).toLocaleString("sv").replace(" ", "_").replaceAll("-", "_").replaceAll(":", "_");
                        e.href = URL.createObjectURL(n), e.download = "".concat(this.$STRINGS.APP_NAME, "_").concat(i, ".json"), 
                        e.click(), URL.revokeObjectURL(e.href);
                    },
                    resetProperties: function() {
                        this.ticket = "", this.error = "", this.mode = "", this.hits = null, this.selectedDatabases = 0, 
                        this.tableMode = 0;
                    },
                    fetchData: function(t) {
                        this.resetProperties(), this.hits = [];
                        var n, e = Tt(t);
                        try {
                            for (e.s(); !(n = e.n()).done; ) {
                                var i = n.value;
                                this.hits.push(this.parseResults(i));
                            }
                        } catch (t) {
                            e.e(t);
                        } finally {
                            e.f();
                        }
                    }
                }
            };
            e(2556), e(8973);
            var Et = e(5893), Dt = e(5255), Ot = e(4506), Rt = (0, M.Z)(Lt, g, [], !1, null, "54679682", null);
            Q()(Rt, {
                VAppBar: ft.Z,
                VAppBarTitle: Ct.Z,
                VBtn: V.Z,
                VCard: Et.Z,
                VCardTitle: Dt.EB,
                VContainer: H.Z,
                VFileInput: Ot.Z,
                VFlex: Y.Z,
                VIcon: W.Z,
                VLayout: J.Z,
                VSpacer: xt.Z,
                VTab: K.Z,
                VTabs: X.Z,
                VToolbarItems: It.lj
            }), Rt.options.__file = "frontend/ResultLocal.vue";
            const kt = {
                components: {
                    ResultLocal: Rt.exports
                }
            };
            var Bt = e(1095), Pt = e(5091), zt = (0, M.Z)(kt, p, [], !1, null, null, null);
            Q()(zt, {
                VApp: Bt.Z,
                VMain: Pt.Z
            }), zt.options.__file = "frontend/AppLocal.vue";
            const Zt = zt.exports;
            e(654);
            i.Z.use(a.Z), i.Z.use(u);
            var _t = {
                mmseqs: e(8615).Z,
                foldseek: e(5473).Z
            };
            window.document.title = _t.foldseek.APP_NAME + " Search Server";
            var Ut = window.matchMedia("(prefers-color-scheme: dark)"), jt = new a.Z({
                icons: {
                    iconfont: "mdiSvg"
                },
                theme: {
                    dark: Ut.matches
                }
            });
            Ut.addEventListener("change", (function(t) {
                jt.framework.theme.dark = t.matches;
            })), i.Z.use({
                install: function(t, n) {
                    t.prototype.$APP = "foldseek", t.prototype.$STRINGS = _t.foldseek, t.prototype.$ELECTRON = !1, 
                    t.prototype.$LOCAL = !0, t.prototype.$MDI = {
                        History: h.BBX,
                        ChevronLeft: h.gAv,
                        ChevronRight: h.zrb,
                        ClockOutline: h.R1X,
                        AlertCircleOutline: h._gM,
                        HelpCircleOutline: h.Gir,
                        Magnify: h.I0v,
                        Tune: h.S3d,
                        Dns: h.cfj,
                        ReorderHorizontal: h.Qjn,
                        Delete: h.x9U,
                        FileDownloadOutline: h.wLz,
                        CloudDownloadOutline: h.REA,
                        FormatListBulleted: h.Ir0,
                        Label: h.KB_,
                        LabelOutline: h.iz_,
                        NotificationClearAll: h.Tal,
                        ProgressWrench: h.Oy8,
                        Restore: h.mBz,
                        Fullscreen: h.h40,
                        ArrowRightCircle: h.BzZ,
                        ArrowRightCircleOutline: h.LHZ,
                        Circle: h.mdD,
                        CircleHalf: h.dMH,
                        PlusBox: h.U1m,
                        MinusBox: h.PeF
                    }, t.prototype.__OS__ = {
                        arch: "web",
                        platform: "web"
                    }, t.prototype.mmseqsVersion = "web", t.prototype.saveResult = function() {}, t.prototype.handleTitleBarDoubleClick = function() {};
                }
            });
            new i.Z({
                el: "#app",
                vuetify: jt,
                render: function(t) {
                    return t(Zt);
                }
            });
        },
        9837: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, 'body, svg text, #app.electron {\n    font-family: system-ui, -apple-system, BlinkMacSystemFont, \'Segoe UI\', Roboto, Oxygen, Ubuntu, Cantarell, \'Open Sans\', \'Helvetica Neue\', sans-serif !important;\n}\n\nbody {\n    background-color: #fff;\n}\n\n@media screen and (prefers-color-scheme: dark) {\n    html, body {\n        background-color: #121212;\n        color-scheme: dark;\n    }\n}\n\nsvg a {\n    cursor: pointer;\n}\n\n.monospace, .mono, pre {\n    font-family: ui-monospace, Inconsolata, Consolas, Menlo, Monaco, "Cascadia Mono", "Segoe UI Mono", "Roboto Mono", "Oxygen Mono", "Ubuntu Monospace", "Source Code Pro", "Fira Mono", "Droid Sans Mono", "Courier New", monospace;\n}\n\n.loading {\n    -webkit-animation: spin 1000ms infinite linear;\n    animation: spin 1000ms infinite linear;\n}\n\n@-webkit-keyframes spin {\n    0% {\n        -webkit-transform: rotate(0deg);\n        transform: rotate(0deg);\n    }\n    100% {\n        -webkit-transform: rotate(359deg);\n        transform: rotate(359deg);\n    }\n}\n@keyframes spin {\n    0% {\n        -webkit-transform: rotate(0deg);\n        transform: rotate(0deg);\n    }\n    100% {\n        -webkit-transform: rotate(359deg);\n        transform: rotate(359deg);\n    }\n}\n\n.input-group .tooltip label {\n    max-width: 100%;\n}\n\nmain.content {\n    max-width: 1536px;\n}\n\n@media print {\n    nav.v-navigation-drawer, header.v-app-bar {\n        display: none !important;\n    }\n    main {\n        padding: 1cm !important;\n    }\n    .v-card, .v-sheet {\n        border: 0px solid transparent !important;\n        outline: 0px solid transparent !important;\n        box-shadow: none !important;\n    }\n}\n\n#app.electron a {\n    -webkit-user-drag: none;\n}\n\n#app.electron .v-toolbar__content, #app.electron .v-input label {\n    user-select: none;\n}', "", {
                version: 3,
                sources: [ "webpack://./frontend/assets/style.css" ],
                names: [],
                mappings: "AAAA;IACI,8JAA8J;AAClK;;AAEA;IACI,sBAAsB;AAC1B;;AAEA;IACI;QACI,yBAAyB;QACzB,kBAAkB;IACtB;AACJ;;AAEA;IACI,eAAe;AACnB;;AAEA;IACI,gOAAgO;AACpO;;AAEA;IACI,8CAA8C;IAC9C,sCAAsC;AAC1C;;AAEA;IACI;QACI,+BAA+B;QAC/B,uBAAuB;IAC3B;IACA;QACI,iCAAiC;QACjC,yBAAyB;IAC7B;AACJ;AACA;IACI;QACI,+BAA+B;QAC/B,uBAAuB;IAC3B;IACA;QACI,iCAAiC;QACjC,yBAAyB;IAC7B;AACJ;;AAEA;IACI,eAAe;AACnB;;AAEA;IACI,iBAAiB;AACrB;;AAEA;IACI;QACI,wBAAwB;IAC5B;IACA;QACI,uBAAuB;IAC3B;IACA;QACI,wCAAwC;QACxC,yCAAyC;QACzC,2BAA2B;IAC/B;AACJ;;AAEA;IACI,uBAAuB;AAC3B;;AAEA;IACI,iBAAiB;AACrB",
                sourcesContent: [ 'body, svg text, #app.electron {\n    font-family: system-ui, -apple-system, BlinkMacSystemFont, \'Segoe UI\', Roboto, Oxygen, Ubuntu, Cantarell, \'Open Sans\', \'Helvetica Neue\', sans-serif !important;\n}\n\nbody {\n    background-color: #fff;\n}\n\n@media screen and (prefers-color-scheme: dark) {\n    html, body {\n        background-color: #121212;\n        color-scheme: dark;\n    }\n}\n\nsvg a {\n    cursor: pointer;\n}\n\n.monospace, .mono, pre {\n    font-family: ui-monospace, Inconsolata, Consolas, Menlo, Monaco, "Cascadia Mono", "Segoe UI Mono", "Roboto Mono", "Oxygen Mono", "Ubuntu Monospace", "Source Code Pro", "Fira Mono", "Droid Sans Mono", "Courier New", monospace;\n}\n\n.loading {\n    -webkit-animation: spin 1000ms infinite linear;\n    animation: spin 1000ms infinite linear;\n}\n\n@-webkit-keyframes spin {\n    0% {\n        -webkit-transform: rotate(0deg);\n        transform: rotate(0deg);\n    }\n    100% {\n        -webkit-transform: rotate(359deg);\n        transform: rotate(359deg);\n    }\n}\n@keyframes spin {\n    0% {\n        -webkit-transform: rotate(0deg);\n        transform: rotate(0deg);\n    }\n    100% {\n        -webkit-transform: rotate(359deg);\n        transform: rotate(359deg);\n    }\n}\n\n.input-group .tooltip label {\n    max-width: 100%;\n}\n\nmain.content {\n    max-width: 1536px;\n}\n\n@media print {\n    nav.v-navigation-drawer, header.v-app-bar {\n        display: none !important;\n    }\n    main {\n        padding: 1cm !important;\n    }\n    .v-card, .v-sheet {\n        border: 0px solid transparent !important;\n        outline: 0px solid transparent !important;\n        box-shadow: none !important;\n    }\n}\n\n#app.electron a {\n    -webkit-user-drag: none;\n}\n\n#app.electron .v-toolbar__content, #app.electron .v-input label {\n    user-select: none;\n}' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        5426: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, '\n.residues {\n    font-family: InconsolataClustal, Inconsolata, Consolas, Menlo, Monaco, "Cascadia Mono", "Segoe UI Mono", "Roboto Mono", "Oxygen Mono", "Ubuntu Monospace", "Source Code Pro", "Fira Mono", "Droid Sans Mono", "Courier New", monospace;\n    white-space: pre;\n}\n.alignment-wrapper-inner {\n    display: inline-block;\n    overflow-x: auto;\n}\n.alignment-wrapper-inner .line {\n    display: inline-block;\n    margin-bottom: 0.5em;\n    white-space: nowrap;\n}\n', "", {
                version: 3,
                sources: [ "webpack://./frontend/Alignment.vue" ],
                names: [],
                mappings: ";AA6FA;IACA,sOAAA;IACA,gBAAA;AACA;AACA;IACA,qBAAA;IACA,gBAAA;AACA;AACA;IACA,qBAAA;IACA,oBAAA;IACA,mBAAA;AACA",
                sourcesContent: [ '<template>\n    <div class="alignment-wrapper-inner">\n        <span class="monospace" v-for="i in Math.max(1, Math.ceil(alignment.alnLength / lineLen))" :key="i">\n            <span class="line">\n                Q&nbsp;{{padNumber(getQueryRowStartPos(i), (Math.max(alignment.qStartPos, alignment.dbStartPos) + alignment.alnLength+"").length, \'&nbsp;\')}}&nbsp;<span class="residues">{{alignment.qAln.substring((i-1)*lineLen,  (i-1)*lineLen+lineLen)}}</span>\n                <br>\n                {{\'&nbsp;\'.repeat(3+(Math.max(alignment.qStartPos, alignment.dbStartPos) + alignment.alnLength+"").length)}}<span class="residues">{{formatAlnDiff(alignment.qAln.substring((i-1)*lineLen,  (i-1)*lineLen+lineLen), alignment.dbAln.substring((i-1)*lineLen, (i-1)*lineLen+lineLen))}}</span>\n                <br>\n                T&nbsp;{{padNumber(getTargetRowStartPos(i), (Math.max(alignment.qStartPos, alignment.dbStartPos) + alignment.alnLength+"").length, \'&nbsp;\')}}&nbsp;<span class="residues" @pointerup="onSelectText(i)">{{alignment.dbAln.substring((i-1)*lineLen, (i-1)*lineLen+lineLen)}}</span>\n            </span><br>\n        </span>\n        <small v-if="$APP == \'foldseek\'" style="float:right">Select target residues to highlight their structure</small>\n    </div>\n</template>\n\n<script>\n\n// cat blosum62.out  | grep -v \'^#\' | awk \'NR == 1 { for (i = 1; i <= NF; i++) { r[i] = $i; } next; } { col = $1; for (i = 2; i <= NF; i++) { print col,r[i-1],$i; } }\' | awk \'$3 > 0 && $1 != $2 { printf "\\""$1""$2"\\",";}\'\nconst blosum62Sim = [\n    "AG", "AS", "DE", "DN",\n    "ED", "EK", "EQ", "FL",\n    "FM", "FW", "FY", "GA",\n    "HN", "HQ", "HY", "IL",\n    "IM", "IV", "KE", "KQ",\n    "KR", "LF", "LI", "LM",\n    "LV", "MF", "MI", "ML",\n    "MV", "ND", "NH", "NQ",\n    "NS", "QE", "QH", "QK",\n    "QN", "QR", "RK", "RQ",\n    "SA", "SN", "ST", "TS",\n    "VI", "VL", "VM", "WF",\n    "WY", "YF", "YH", "YW"\n]\n\n// Get the first and last non-null values in a map between a range\nfunction getRange(map, start, end) {\n    let first = null, last = null\n    for (let i = start; i <= end; i++) {\n\tlet val = map[i]\n\tif (val !== null) {\n\t    if (first === null) first = val\n\t    last = val\n\t}\n    }\n    return [first, last]\n}\n\nexport default {\n    props: [\'alignment\', \'lineLen\', \'queryMap\', \'targetMap\'],\n    methods: {\n        // Get the index of a given residue in the alignment\n        getQueryIndex(index) { return this.queryMap[index] },\n        getTargetIndex(index) { return this.targetMap[index] },\n        getFirstResidueNumber(map, i) {\n            let start = this.lineLen * (i - 1)\n            while (map[start] === null) start--\n            return map[start]\n        },\n        getQueryRowStartPos(i) { return this.getFirstResidueNumber(this.queryMap, i) },\n        getTargetRowStartPos(i) { return this.getFirstResidueNumber(this.targetMap, i) },\n        formatAlnDiff(seq1, seq2) {\n            if (seq1.length != seq2.length) return \'\'\n            var res = \'\'\n            for (var i = 0; i < seq1.length; i++) {\n                if (seq1[i] == seq2[i]) res += seq1[i];\n                else if (blosum62Sim.indexOf(seq1[i] + seq2[i]) != -1) res += \'+\';\n                else res += \' \';\n            }\n            return res;\n        },\n        padNumber(nr, n, str){\n            return Array(n - String(nr).length + 1).join(str || \'0\') + nr\n        },\n        onSelectText(i) {\n            var selection = window.getSelection()\n\n            // In case of backwards selection\n            var [offsetStart, offsetEnd] = [\n                selection.anchorOffset, selection.focusOffset\n            ].sort((a, b) => a - b)\n\n            var length = offsetEnd - offsetStart\n            var relStart = (i - 1) * this.lineLen + offsetStart\n            var relEnd = relStart + length - 1 // the selection is inclusive\n\n            var [start, end] = getRange(this.targetMap, relStart, relEnd)\n            this.$emit(\'selected\', [start, end])\n        }\n    }, \n}\n<\/script>\n\n<style>\n.residues {\n    font-family: InconsolataClustal, Inconsolata, Consolas, Menlo, Monaco, "Cascadia Mono", "Segoe UI Mono", "Roboto Mono", "Oxygen Mono", "Ubuntu Monospace", "Source Code Pro", "Fira Mono", "Droid Sans Mono", "Courier New", monospace;\n    white-space: pre;\n}\n.alignment-wrapper-inner {\n    display: inline-block;\n    overflow-x: auto;\n}\n.alignment-wrapper-inner .line {\n    display: inline-block;\n    margin-bottom: 0.5em;\n    white-space: nowrap;\n}\n</style>\n' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        6696: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, "\n.alignment-wrapper-outer {\n    display: inline-flex;\n    flex-direction: row;\n    flex-wrap: nowrap;\n    justify-content: center;\n    align-items: stretch;\n    width: 100%;\n}\n.alignment-wrapper-inner {\n    flex: 2;\n    margin: auto;\n    display: flex;\n    flex-direction: column;\n    align-items: end;\n}\n.alignment-structure-wrapper {\n    flex: 1;\n    min-width:450px;\n    margin: 0;\n    margin-bottom: auto;\n}\n@media screen and (max-width: 960px) {\n.alignment-wrapper-outer {\n        display: flex;\n        flex-direction: column;\n}\n.alignment-structure-wrapper {\n        padding-top: 1em;\n}\n}\n@media screen and (min-width: 961px) {\n.alignment-structure-wrapper {\n        padding-left: 2em;\n}\n}\n\n", "", {
                version: 3,
                sources: [ "webpack://./frontend/AlignmentPanel.vue" ],
                names: [],
                mappings: ";AA2EA;IACA,oBAAA;IACA,mBAAA;IACA,iBAAA;IACA,uBAAA;IACA,oBAAA;IACA,WAAA;AACA;AACA;IACA,OAAA;IACA,YAAA;IACA,aAAA;IACA,sBAAA;IACA,gBAAA;AACA;AAEA;IACA,OAAA;IACA,eAAA;IACA,SAAA;IACA,mBAAA;AACA;AAEA;AACA;QACA,aAAA;QACA,sBAAA;AACA;AACA;QACA,gBAAA;AACA;AACA;AAEA;AACA;QACA,iBAAA;AACA;AACA",
                sourcesContent: [ '<template>\n    <div class="alignment-wrapper-outer" slot="content">\n        <Alignment\n            :key="`aln2-${alignment.id}`"\n            :alignment="alignment"\n            :lineLen="lineLen"\n            :queryMap="queryMap"\n            :targetMap="targetMap"\n            @selected="setUserSelection"\n        />\n        <div v-if="$APP == \'foldseek\'" class="alignment-structure-wrapper">\n            <StructureViewer\n                :key="`struc2-${alignment.id}`"\n                :alignment="alignment"\n                :queryMap="queryMap"\n                :targetMap="targetMap"\n                :hits="hits"\n                bgColorLight="white"\n                bgColorDark="#1E1E1E"\n                qColor="lightgrey"\n                tColor="red"\n                qRepr="cartoon"\n                tRepr="cartoon"\n                ref="structureViewer"\n            />\n        </div>\n    </div>\n</template>\n\n<script>\nimport Alignment from \'./Alignment.vue\'\n\n// Map 0-based indices in the alignment to corresponding 1-based indices in the structure\nfunction makePositionMap(realStart, alnString) {\n    let map = Array(alnString.length);\n    for (let i = 0, gaps = 0; i < alnString.length; i++) {\n        if (alnString[i] === \'-\') {\n            map[i] = null;\n            gaps++;\n        } else {\n            map[i] = realStart + i - gaps;\n        }\n    }\n    return map\n}\n\nexport default {\n    components: { StructureViewer: () => __APP__ == "foldseek" ? import(\'./StructureViewer.vue\') : null, Alignment },\n    data: () => ({\n        queryMap: null,\n        targetMap: null,\n    }),\n    props: {\n        alignment: { type: Object, required: true, },\n        lineLen: { type: Number, required: true, },\n        hits: { type: Object }\n    },\n    methods: {\n        setUserSelection([start, end]) {\n            if (!this.alignment) return\n            if (__APP__ != "foldseek") return\n            this.$refs.structureViewer.setSelectionData(start, end)\n        },\n        updateMaps() {\n            if (!this.alignment) return\n            this.queryMap = makePositionMap(this.alignment.qStartPos, this.alignment.qAln)\n            this.targetMap = makePositionMap(this.alignment.dbStartPos, this.alignment.dbAln)\n        },\n    },\n    watch: { \'alignment\': function() { this.updateMaps() } },\n    beforeMount() { this.updateMaps() },\n}\n<\/script>\n\n<style>\n.alignment-wrapper-outer {\n    display: inline-flex;\n    flex-direction: row;\n    flex-wrap: nowrap;\n    justify-content: center;\n    align-items: stretch;\n    width: 100%;\n}\n.alignment-wrapper-inner {\n    flex: 2;\n    margin: auto;\n    display: flex;\n    flex-direction: column;\n    align-items: end;\n}\n\n.alignment-structure-wrapper {\n    flex: 1;\n    min-width:450px;\n    margin: 0;\n    margin-bottom: auto;\n}\n\n@media screen and (max-width: 960px) {\n    .alignment-wrapper-outer {\n        display: flex;\n        flex-direction: column;\n    }\n    .alignment-structure-wrapper {\n        padding-top: 1em;\n    }\n}\n\n@media screen and (min-width: 961px) {\n    .alignment-structure-wrapper {\n        padding-left: 2em;\n    }\n}\n\n</style>\n' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        8260: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, "\n[data-v-5976e89a] .v-app-bar-title__content {\n    text-overflow: revert !important;\n}\n[data-v-5976e89a] .theme--light.v-navigation-drawer {\n    background-color: #f5f5f5;\n    border-color: #f5f5f5;\n    /* transition-duration: 0s !important; */\n    /* transition-timing-function: linear; */\n}\n[data-v-5976e89a] .theme--dark.v-navigation-drawer {\n    background-color: #212121;\n    border-color: #212121;\n}\n", "", {
                version: 3,
                sources: [ "webpack://./frontend/Navigation.vue" ],
                names: [],
                mappings: ";AA+FA;IACA,gCAAA;AACA;AACA;IACA,yBAAA;IACA,qBAAA;IACA,wCAAA;IACA,wCAAA;AACA;AAEA;IACA,yBAAA;IACA,qBAAA;AACA",
                sourcesContent: [ '<template>\n<div>\n<v-navigation-drawer v-if="!$LOCAL" stateless app permanent clipped :mini-variant="mini" :expand-on-hover="false" ref="drawer">\n    <v-list v-if="!$LOCAL">\n        <v-list-item to="/search">\n            <v-list-item-action>\n                <v-icon>{{ $MDI.Magnify }}</v-icon>\n            </v-list-item-action>\n            <v-list-item-content>\n                <v-list-item-title>Search</v-list-item-title>\n            </v-list-item-content>\n        </v-list-item>\n        \n        <v-list-item>\n            <v-list-item-action>\n                <v-icon>{{ $MDI.Upload }}</v-icon>\n            </v-list-item-action>\n            <v-list-item-content>\n                <v-list-item-title>Upload JSON</v-list-item-title>\n            </v-list-item-content>\n        </v-list-item>\n\n        <router-view name="sidebar"></router-view>\n        <history v-if="!$LOCAL" />\n\n        <v-list-item v-if="$ELECTRON" to="/preferences">\n            <v-list-item-action>\n                <v-icon>{{ $MDI.Tune }}</v-icon>\n            </v-list-item-action>\n            <v-list-item-content>\n                <v-list-item-title>Preferences</v-list-item-title>\n            </v-list-item-content>\n        </v-list-item>\n    </v-list>\n</v-navigation-drawer>\n<v-app-bar v-on:dblclick.native="electronHandleTitleBarDoubleClick()" app :height="$ELECTRON ? \'72px\' : \'48px\'" fixed clipped-left :class="[\'ml-0\', \'pl-3\', $ELECTRON ? \'pt-2\' : null]" :style="{\'-webkit-app-region\': $ELECTRON ? \'drag\' : null, \'-webkit-user-select\': $ELECTRON ? \'none\' : null}">\n    <v-app-bar-nav-icon v-if="!$LOCAL" :input-value="!mini ? \'activated\' : undefined" @click.stop="toggleMini"></v-app-bar-nav-icon>\n    <v-app-bar-title>\n        <router-link v-if="!$LOCAL" to="/" style="color: inherit; text-decoration: none">{{ $STRINGS.APP_NAME }} Search</router-link>\n        <span v-if="$LOCAL">{{ $STRINGS.APP_NAME }} Search</span>\n    </v-app-bar-title>\n    <object style="margin-left:8px; display: inline-block; width: 38px;height: 38px;vertical-align: middle"\n            v-if="$APP == \'mmseqs\'"\n            type="image/svg+xml"\n            data="./assets/marv1.svg"\n            aria-hidden="true">\n        <img src="./assets/marv1.png" style="max-width:100%" />\n    </object>\n    <img v-if="$APP == \'foldseek\'" src="./assets/marv-foldseek-small.png" style="margin-left:8px; display: inline-block; width: 48px;height: 48px;vertical-align: middle" aria-hidden="true" />\n\n    <v-spacer></v-spacer>\n    <v-toolbar-items v-once v-if="!$ELECTRON" class="hidden-sm-and-down">\n        <v-btn text rel="external noopener" target="_blank"\n               v-for="i in ($STRINGS.NAV_URL_COUNT - 0)" :key="i" :href="$STRINGS[\'NAV_URL_\' + i]">{{ $STRINGS["NAV_TITLE_" + i]}}</v-btn>\n    </v-toolbar-items>\n</v-app-bar>\n\n</div>\n</template>\n\n<script>\nimport History from \'./History.vue\';\n\nexport default {\n    components : { History },\n    data: () => ({\n        mini: true,\n    }),\n    created() {\n        this.$root.$on(\'multi\', this.shouldExpand);\n    },\n    mounted() {\n        // defeat https://github.com/vuetifyjs/vuetify/pull/14523\n        if (!__LOCAL__) Object.defineProperty(this.$refs.drawer._data, \'isMouseover\', { get: () => { false } });\n        \n        console.log(this)\n    },\n    beforeDestroy() {\n        this.$root.$off(\'multi\', this.shouldExpand);\n    },\n    methods: {\n        shouldExpand(expand) {\n            this.mini = !expand;\n        },\n        toggleMini() {\n            this.mini = !this.mini;\n        },\n        electronHandleTitleBarDoubleClick() {\n            this.handleTitleBarDoubleClick();\n        }\n    }\n}\n<\/script>\n\n<style scoped>\n::v-deep .v-app-bar-title__content {\n    text-overflow: revert !important;\n}\n::v-deep .theme--light.v-navigation-drawer {\n    background-color: #f5f5f5;\n    border-color: #f5f5f5;\n    /* transition-duration: 0s !important; */\n    /* transition-timing-function: linear; */\n}\n\n::v-deep .theme--dark.v-navigation-drawer {\n    background-color: #212121;\n    border-color: #212121;\n}\n</style>\n' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        4569: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => p
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r), o = e(1667), l = e.n(o), A = new URL(e(42), e.b), c = new URL(e(901), e.b), d = s()(a()), u = l()(A), h = l()(c);
            d.push([ t.id, "\n.panel-root[data-v-0d9b5935], .panel-content[data-v-0d9b5935] {\n    flex-direction: column;\n}\n.panel-root header[data-v-0d9b5935], .panel-content[data-v-0d9b5935] {\n    contain: content;\n}\n.panel-root nav[data-v-0d9b5935] {\n    flex: 0;\n}\n.panel-root .force-fill-height[data-v-0d9b5935] {\n    display: flex;\n    height: 100% !important;\n}\n.panel-root[data-v-0d9b5935] .v-toolbar {\n    background-repeat: repeat;\n}\n.theme--light .panel-root[data-v-0d9b5935] .v-toolbar {\n    background: url(" + u + ");\n}\n.theme--dark .panel-root[data-v-0d9b5935] .v-toolbar {\n    background: url(" + h + ");\n}\n.panel-root[data-v-0d9b5935] .text-h6 {\n    margin-bottom: -5px;\n}\n.panel-root[data-v-0d9b5935] .text-h6 i.v-icon {\n    font-size: 1em;\n    vertical-align: bottom;\n}\n", "", {
                version: 3,
                sources: [ "webpack://./frontend/Panel.vue" ],
                names: [],
                mappings: ";AAsDA;IACA,sBAAA;AACA;AAEA;IACA,gBAAA;AACA;AAEA;IACA,OAAA;AACA;AAEA;IACA,aAAA;IACA,uBAAA;AACA;AAEA;IACA,yBAAA;AACA;AAEA;IACA,mDAAA;AAEA;AAEA;IACA,mDAAA;AACA;AAEA;IACA,mBAAA;AACA;AAEA;IACA,cAAA;IACA,sBAAA;AACA",
                sourcesContent: [ "<template>\n    <div :class=\"['panel-root', elevation != null ? 'elevation-' + elevation : null ]\">\n        <v-toolbar v-if=\"!!$slots['header'] || !!header\" text dense dark>\n            <v-btn v-if=\"collapsible\" style=\"margin-top:0;margin-left:-15px;\" icon plain  @click=\"isCollapsed = !isCollapsed\" :aria-expanded=\"isCollapsed ? 'false' : 'true'\" :aria-controls=\"uuid\">\n                <v-icon v-if=\"isCollapsed\">\n                    {{ $MDI.PlusBox }}\n                </v-icon>\n                <v-icon v-else>\n                    {{ $MDI.MinusBox }}\n                </v-icon>\n            </v-btn>\n            <span class=\"text-h6 align-end\">\n                <slot v-if=\"$slots['header']\" name=\"header\"></slot>\n                <template v-else>{{ header }}</template>\n            </span>\n            <v-spacer></v-spacer>\n            <slot name=\"toolbar-extra\"></slot>\n        </v-toolbar>\n        <v-card rounded=\"0\" :class=\"['panel', { 'd-flex' : flex }, { 'force-fill-height' : fillHeight }]\" v-if=\"!isCollapsed\" :id=\"uuid\">\n            <v-card-text v-if=\"$slots['desc']\" class=\"subheading justify\">\n                <slot name=\"desc\"></slot>\n            </v-card-text>\n            <v-card-text v-if=\"$slots['content']\" :class=\"['panel-content', 'justify', { 'd-flex' : flex }]\">\n                <slot name=\"content\"></slot>\n            </v-card-text>\n        </v-card>\n    </div>\n</template>\n\n<script>\nlet uuid = 0;\nexport default {\n    name: 'panel',\n    props: { \n        header : { default: '', type: String }, \n        'fillHeight' : { default: false, type: Boolean }, \n        'collapsible' : { default: false, type: Boolean },\n        'collapsed' : { default: false, type: Boolean },\n        'flex' : { default: true, type: Boolean },\n        'elevation' : { default: null, type: Number }\n    },\n    data() {\n        return {\n            isCollapsed: this.collapsed,\n        }\n    },\n    beforeCreate() {\n        this.uuid = 'panel-' + uuid.toString();\n        uuid += 1;\n    },\n}\n<\/script>\n\n<style scoped>\n.panel-root, .panel-content {\n    flex-direction: column;\n}\n\n.panel-root header, .panel-content {\n    contain: content;\n}\n\n.panel-root nav {\n    flex: 0;\n}\n\n.panel-root .force-fill-height {\n    display: flex;\n    height: 100% !important;\n}\n\n.panel-root >>> .v-toolbar {\n    background-repeat: repeat;\n}\n\n.theme--light .panel-root >>> .v-toolbar {\n    background: url('./assets/spiration-dark.png');\n    \n}\n\n.theme--dark .panel-root >>> .v-toolbar {\n    background: url('./assets/spiration-darker.png');\n}\n\n.panel-root >>> .text-h6 {\n    margin-bottom: -5px;\n}\n\n.panel-root >>> .text-h6 i.v-icon {\n    font-size: 1em;\n    vertical-align: bottom;\n}\n</style>" ],
                sourceRoot: ""
            } ]);
            const p = d;
        },
        864: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, "\n[data-v-54679682] .v-app-bar-title__content {\n    text-overflow: revert !important;\n}\n", "", {
                version: 3,
                sources: [ "webpack://./frontend/ResultLocal.vue" ],
                names: [],
                mappings: ";AAoKA;IACA,gCAAA;AACA",
                sourcesContent: [ '<template>\n    <div>\n        <v-app-bar app :height="\'48px\'" fixed clipped-left>\n            <img height="28px" src="data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbDpzcGFjZT0icHJlc2VydmUiIHN0eWxlPSJmaWxsLXJ1bGU6ZXZlbm9kZDtjbGlwLXJ1bGU6ZXZlbm9kZDtzdHJva2UtbGluZWNhcDpyb3VuZDtzdHJva2UtbGluZWpvaW46cm91bmQ7c3Ryb2tlLW1pdGVybGltaXQ6MTAiIHZpZXdCb3g9IjAgMCA0NjggMzA2Ij48cGF0aCBkPSJNMzcyIDIwMnMxNC0xIDM3LTE5YzIzLTE3IDQwLTQ5IDU1LTU1bC0xMTQgMjQtNCAzMiAyNiAxOFoiIHN0eWxlPSJmaWxsOiNmN2QxOGE7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik02MiAxMzlTODcgMjEgMjY5IDJsMSAxLTQ2IDYxcy00MC0zLTU1IDdjMCAwIDE5LTEzIDY5LTRzNTAtMjAgNTAtMjAgOCAyMiAwIDI5bDI5IDE0LTE4IDRzMTI1LTEyIDE2NyAzM2MwIDAtMjYgMTctNjAgMjAtNTYgNS02MiAyMi02MiAyMnMyNS0xMCA0MyA0bC0yMiA5czE1IDggMTUgMjNsLTI2IDEwczM2LTE4IDUyLTdsLTI0IDE4czIzIDMgMzggMTVsLTMyIDhzMTUgMiAyNyAzMWwtNDUtNnM3IDkgNCAzMGwtMjUtMjJzLTE3IDQ2LTE1OCAyQzQ5IDI0MCA1NiAyMjEgNTAgMTkxbC0yNi0xczItMTUgMTgtMjFMMiAxNDJzMjQtMTMgNDItOGwtOC0yNXMyOSAxMSAyNiAzMFoiIHN0eWxlPSJmaWxsOiNlMTMyMTM7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0xMDEgMjUzYy00Ni0yMyA4LTEzNCAzNy0xNTEgMjgtMTYgNTcgNyA2MyAxOSAwIDAgMjMtMTggNTctN3M0OSA0NyAzNiAxMTVjLTggNDEtMjQgNTgtMzUgNjUtNyA0LTE0IDUtMjEgMy0yNS02LTEwNS0yNy0xMzctNDRaIiBzdHlsZT0iZmlsbDojZjdkMThhO2ZpbGwtcnVsZTpub256ZXJvO3N0cm9rZTojMDAwO3N0cm9rZS13aWR0aDo0LjQ4cHgiLz48cGF0aCBkPSJNMTM2IDExMnMtNDEtMTAtNTYgMThjLTE1IDI3IDEyIDM4IDI3IDQzIDE2IDQgNDcgNCA1Ny0xM3MtMS0zOC0yOC00OFoiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0xMTYgMTYwYzE2IDggMzQtMzcgMjAtNDQtMTQtNi00MCAzNS0yMCA0NFoiIHN0eWxlPSJmaWxsLXJ1bGU6bm9uemVybztzdHJva2U6IzAwMDtzdHJva2Utd2lkdGg6NC40OHB4Ii8+PHBhdGggZD0iTTI4NCAxNDhjLTQxLTE1LTU5IDUtNjUgMjJzMiA0NCA0MiA1MyA1MC00IDU2LTE5YzUtMTYgNi00MS0zMy01NloiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0yNDggMTk5YzE5IDkgNDctNDEgMjMtNTJzLTQzIDQzLTIzIDUyWm0tODUtMTVjMS04IDIwLTEgMjAgNSAwIDctOSA4LTEyIDctNC0xLTktNi04LTEyWiIgc3R5bGU9ImZpbGwtcnVsZTpub256ZXJvO3N0cm9rZTojMDAwO3N0cm9rZS13aWR0aDo0LjQ4cHgiLz48cGF0aCBkPSJNMTMyIDEyMGM3IDMtMiAxNS02IDEyczMtMTQgNi0xMlptMTI4IDMwYzcgMy0yIDE1LTYgMTItNC0yIDMtMTQgNi0xMloiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm8iLz48cGF0aCBkPSJtMTE1IDIxMiA5LTRzLTggNyAwIDEzYzggNyAyNS00IDQ2LTEgMjEgNCA0MCAxOSA1NSAyMSAxNiAzIDI0IDEgMjMtNC0xLTYgNSA3IDUgNyIgc3R5bGU9ImZpbGw6bm9uZTtmaWxsLXJ1bGU6bm9uemVybztzdHJva2U6IzAwMDtzdHJva2Utd2lkdGg6NC40OHB4Ii8+PC9zdmc+" />\n            &nbsp;\n            <v-app-bar-title class="ml-2">{{ $STRINGS.APP_NAME }} Search</v-app-bar-title>\n            <v-spacer />\n            <v-file-input\n                id="uploadData"\n                class="shrink"\n                type="file"\n                accept="application/json"\n                placeholder="Load JSON data file"\n                style="position: relative; top: 30%;"\n                @change="uploadData" \n                single-line\n                outlined\n                filled\n                flat\n                dense\n            />\n            <v-toolbar-items>\n                <v-btn text @click="downloadData">\n                    <v-icon>\n                        {{ $MDI.FileDownloadOutline }}\n                    </v-icon>\n                </v-btn>\n                <v-btn text rel="external noopener" target="_blank" class="hidden-sm-and-down"\n                       v-for="i in ($STRINGS.NAV_URL_COUNT - 0)" :key="i" :href="$STRINGS[\'NAV_URL_\' + i]">{{ $STRINGS["NAV_TITLE_" + i]}}</v-btn>\n            </v-toolbar-items>\n        </v-app-bar>\n        <v-tabs v-if="hits" center-active grow style="margin-bottom: 1em" show-arrows>\n            <v-tab v-for="(entry, index) in hits" :key="entry.query.header" @click="changeResult(index)">\n                {{ entry.query.header }} ({{ entry.results[0].alignments ? entry.results[0].alignments.length : 0 }})\n            </v-tab>\n        </v-tabs>\n        <ResultView\n            v-if="hits"\n            :key="currentIndex"\n            :ticket="ticket"\n            :error="error"\n            :mode="mode"\n            :hits="currentResult"\n            :selectedDatabases="selectedDatabases"\n            :tableMode="tableMode"\n        />\n        <v-container grid-list-md fluid pa-2 v-else>\n            <v-layout wrap>\n                <v-flex xs12>\n                    <v-card rounded="0">\n                        <v-card-title primary-title class="mb-0 pa-4">\n                            No data loaded\n                        </v-card-title>\n                    </v-card>\n                </v-flex>\n            </v-layout>\n        </v-container>\n        <v-container grid-list-md fluid pa-2>\n            <v-layout wrap>\n                <v-flex xs12>\n                    <v-card rounded="0">\n                    <v-card-title primary-title class="pb-0 mb-0">\n                        <div class="text-h5 mb-0">Reference</div>\n                    </v-card-title>\n                    <v-card-title primary-title class="pt-0 mt-0">\n                        <p class="text-subtitle-2 mb-0" v-html="$STRINGS.CITATION"></p>\n                    </v-card-title>\n                    </v-card>\n                </v-flex>\n            </v-layout>\n        </v-container>\n    </div>\n</template>\n\n<script>\nimport ResultMixin from \'./ResultMixin.vue\';\nimport ResultView from \'./ResultView.vue\';\nimport Navigation from \'./Navigation.vue\';\n\nexport default {\n    name: \'result\',\n    mixins: [ResultMixin],\n    components: { ResultView, Navigation },\n    data() {\n        return {\n            currentIndex: 0\n        };\n    },\n    mounted() {\n        document.onreadystatechange = () => {\n            if (document.readyState == "complete") {\n                let div = document.getElementById("data");\n                if (!div) {\n                    return null;\n                }\n                let data = JSON.parse(div.textContent);\n                this.fetchData(data);\n            }\n        }\n    },\n    computed: {\n        currentResult() {\n            if (this.hits === null)\n                return null;\n            return this.hits[this.currentIndex];\n        },\n        currentQuery() {\n            if (this.hits === null)\n                return "";\n            return this.hits[this.currentIndex].query.header;\n        }\n    },\n    methods: {\n        changeResult(newRes) {\n            this.currentIndex = newRes;\n            this.setColorScheme();\n        },\n        uploadData(file) {\n            if (!file) {\n                return;\n            }\n            let fr = new FileReader();\n            fr.addEventListener(\n                "load",\n                (e) => {\n                    let data = JSON.parse(e.target.result);\n                    this.fetchData(data);\n                }\n            );\n            fr.readAsText(file)\n        },\n        downloadData() {\n            if (!this.hits) {\n                return null;\n            }\n            const json = JSON.stringify(this.hits);\n            const blob = new Blob([json], { type: \'application/json\' });\n            const link = document.createElement(\'a\');\n            const date = new Date().toLocaleString(\'sv\').replace(\' \', \'_\').replaceAll(\'-\', \'_\').replaceAll(\':\', \'_\');\n            link.href = URL.createObjectURL(blob);\n            link.download = `${this.$STRINGS.APP_NAME}_${date}.json`;\n            link.click();\n            URL.revokeObjectURL(link.href);\n        },\n        resetProperties() {\n            this.ticket = "";\n            this.error = "";\n            this.mode = "";\n            this.hits = null;\n            this.selectedDatabases = 0;\n            this.tableMode = 0;\n        },\n        fetchData(data) {\n            this.resetProperties();\n            this.hits = [];\n            for (let result of data) {\n                this.hits.push(this.parseResults(result));\n            }\n        }\n    }\n};\n<\/script>\n\n<style scoped>\n::v-deep .v-app-bar-title__content {\n    text-overflow: revert !important;\n}\n</style>\n<style>\n.theme--light .panel-root .v-toolbar {\n    background-color: #454545 !important;\n}\n\n.theme--dark .panel-root .v-toolbar {\n    background-color: #1e1e1e !important;\n}\n</style>' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        8742: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, "\n.theme--light .panel-root .v-toolbar {\n    background-color: #454545 !important;\n}\n.theme--dark .panel-root .v-toolbar {\n    background-color: #1e1e1e !important;\n}\n", "", {
                version: 3,
                sources: [ "webpack://./frontend/ResultLocal.vue" ],
                names: [],
                mappings: ";AAyKA;IACA,oCAAA;AACA;AAEA;IACA,oCAAA;AACA",
                sourcesContent: [ '<template>\n    <div>\n        <v-app-bar app :height="\'48px\'" fixed clipped-left>\n            <img height="28px" src="data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbDpzcGFjZT0icHJlc2VydmUiIHN0eWxlPSJmaWxsLXJ1bGU6ZXZlbm9kZDtjbGlwLXJ1bGU6ZXZlbm9kZDtzdHJva2UtbGluZWNhcDpyb3VuZDtzdHJva2UtbGluZWpvaW46cm91bmQ7c3Ryb2tlLW1pdGVybGltaXQ6MTAiIHZpZXdCb3g9IjAgMCA0NjggMzA2Ij48cGF0aCBkPSJNMzcyIDIwMnMxNC0xIDM3LTE5YzIzLTE3IDQwLTQ5IDU1LTU1bC0xMTQgMjQtNCAzMiAyNiAxOFoiIHN0eWxlPSJmaWxsOiNmN2QxOGE7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik02MiAxMzlTODcgMjEgMjY5IDJsMSAxLTQ2IDYxcy00MC0zLTU1IDdjMCAwIDE5LTEzIDY5LTRzNTAtMjAgNTAtMjAgOCAyMiAwIDI5bDI5IDE0LTE4IDRzMTI1LTEyIDE2NyAzM2MwIDAtMjYgMTctNjAgMjAtNTYgNS02MiAyMi02MiAyMnMyNS0xMCA0MyA0bC0yMiA5czE1IDggMTUgMjNsLTI2IDEwczM2LTE4IDUyLTdsLTI0IDE4czIzIDMgMzggMTVsLTMyIDhzMTUgMiAyNyAzMWwtNDUtNnM3IDkgNCAzMGwtMjUtMjJzLTE3IDQ2LTE1OCAyQzQ5IDI0MCA1NiAyMjEgNTAgMTkxbC0yNi0xczItMTUgMTgtMjFMMiAxNDJzMjQtMTMgNDItOGwtOC0yNXMyOSAxMSAyNiAzMFoiIHN0eWxlPSJmaWxsOiNlMTMyMTM7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0xMDEgMjUzYy00Ni0yMyA4LTEzNCAzNy0xNTEgMjgtMTYgNTcgNyA2MyAxOSAwIDAgMjMtMTggNTctN3M0OSA0NyAzNiAxMTVjLTggNDEtMjQgNTgtMzUgNjUtNyA0LTE0IDUtMjEgMy0yNS02LTEwNS0yNy0xMzctNDRaIiBzdHlsZT0iZmlsbDojZjdkMThhO2ZpbGwtcnVsZTpub256ZXJvO3N0cm9rZTojMDAwO3N0cm9rZS13aWR0aDo0LjQ4cHgiLz48cGF0aCBkPSJNMTM2IDExMnMtNDEtMTAtNTYgMThjLTE1IDI3IDEyIDM4IDI3IDQzIDE2IDQgNDcgNCA1Ny0xM3MtMS0zOC0yOC00OFoiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0xMTYgMTYwYzE2IDggMzQtMzcgMjAtNDQtMTQtNi00MCAzNS0yMCA0NFoiIHN0eWxlPSJmaWxsLXJ1bGU6bm9uemVybztzdHJva2U6IzAwMDtzdHJva2Utd2lkdGg6NC40OHB4Ii8+PHBhdGggZD0iTTI4NCAxNDhjLTQxLTE1LTU5IDUtNjUgMjJzMiA0NCA0MiA1MyA1MC00IDU2LTE5YzUtMTYgNi00MS0zMy01NloiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm87c3Ryb2tlOiMwMDA7c3Ryb2tlLXdpZHRoOjQuNDhweCIvPjxwYXRoIGQ9Ik0yNDggMTk5YzE5IDkgNDctNDEgMjMtNTJzLTQzIDQzLTIzIDUyWm0tODUtMTVjMS04IDIwLTEgMjAgNSAwIDctOSA4LTEyIDctNC0xLTktNi04LTEyWiIgc3R5bGU9ImZpbGwtcnVsZTpub256ZXJvO3N0cm9rZTojMDAwO3N0cm9rZS13aWR0aDo0LjQ4cHgiLz48cGF0aCBkPSJNMTMyIDEyMGM3IDMtMiAxNS02IDEyczMtMTQgNi0xMlptMTI4IDMwYzcgMy0yIDE1LTYgMTItNC0yIDMtMTQgNi0xMloiIHN0eWxlPSJmaWxsOiNmZmY7ZmlsbC1ydWxlOm5vbnplcm8iLz48cGF0aCBkPSJtMTE1IDIxMiA5LTRzLTggNyAwIDEzYzggNyAyNS00IDQ2LTEgMjEgNCA0MCAxOSA1NSAyMSAxNiAzIDI0IDEgMjMtNC0xLTYgNSA3IDUgNyIgc3R5bGU9ImZpbGw6bm9uZTtmaWxsLXJ1bGU6bm9uemVybztzdHJva2U6IzAwMDtzdHJva2Utd2lkdGg6NC40OHB4Ii8+PC9zdmc+" />\n            &nbsp;\n            <v-app-bar-title class="ml-2">{{ $STRINGS.APP_NAME }} Search</v-app-bar-title>\n            <v-spacer />\n            <v-file-input\n                id="uploadData"\n                class="shrink"\n                type="file"\n                accept="application/json"\n                placeholder="Load JSON data file"\n                style="position: relative; top: 30%;"\n                @change="uploadData" \n                single-line\n                outlined\n                filled\n                flat\n                dense\n            />\n            <v-toolbar-items>\n                <v-btn text @click="downloadData">\n                    <v-icon>\n                        {{ $MDI.FileDownloadOutline }}\n                    </v-icon>\n                </v-btn>\n                <v-btn text rel="external noopener" target="_blank" class="hidden-sm-and-down"\n                       v-for="i in ($STRINGS.NAV_URL_COUNT - 0)" :key="i" :href="$STRINGS[\'NAV_URL_\' + i]">{{ $STRINGS["NAV_TITLE_" + i]}}</v-btn>\n            </v-toolbar-items>\n        </v-app-bar>\n        <v-tabs v-if="hits" center-active grow style="margin-bottom: 1em" show-arrows>\n            <v-tab v-for="(entry, index) in hits" :key="entry.query.header" @click="changeResult(index)">\n                {{ entry.query.header }} ({{ entry.results[0].alignments ? entry.results[0].alignments.length : 0 }})\n            </v-tab>\n        </v-tabs>\n        <ResultView\n            v-if="hits"\n            :key="currentIndex"\n            :ticket="ticket"\n            :error="error"\n            :mode="mode"\n            :hits="currentResult"\n            :selectedDatabases="selectedDatabases"\n            :tableMode="tableMode"\n        />\n        <v-container grid-list-md fluid pa-2 v-else>\n            <v-layout wrap>\n                <v-flex xs12>\n                    <v-card rounded="0">\n                        <v-card-title primary-title class="mb-0 pa-4">\n                            No data loaded\n                        </v-card-title>\n                    </v-card>\n                </v-flex>\n            </v-layout>\n        </v-container>\n        <v-container grid-list-md fluid pa-2>\n            <v-layout wrap>\n                <v-flex xs12>\n                    <v-card rounded="0">\n                    <v-card-title primary-title class="pb-0 mb-0">\n                        <div class="text-h5 mb-0">Reference</div>\n                    </v-card-title>\n                    <v-card-title primary-title class="pt-0 mt-0">\n                        <p class="text-subtitle-2 mb-0" v-html="$STRINGS.CITATION"></p>\n                    </v-card-title>\n                    </v-card>\n                </v-flex>\n            </v-layout>\n        </v-container>\n    </div>\n</template>\n\n<script>\nimport ResultMixin from \'./ResultMixin.vue\';\nimport ResultView from \'./ResultView.vue\';\nimport Navigation from \'./Navigation.vue\';\n\nexport default {\n    name: \'result\',\n    mixins: [ResultMixin],\n    components: { ResultView, Navigation },\n    data() {\n        return {\n            currentIndex: 0\n        };\n    },\n    mounted() {\n        document.onreadystatechange = () => {\n            if (document.readyState == "complete") {\n                let div = document.getElementById("data");\n                if (!div) {\n                    return null;\n                }\n                let data = JSON.parse(div.textContent);\n                this.fetchData(data);\n            }\n        }\n    },\n    computed: {\n        currentResult() {\n            if (this.hits === null)\n                return null;\n            return this.hits[this.currentIndex];\n        },\n        currentQuery() {\n            if (this.hits === null)\n                return "";\n            return this.hits[this.currentIndex].query.header;\n        }\n    },\n    methods: {\n        changeResult(newRes) {\n            this.currentIndex = newRes;\n            this.setColorScheme();\n        },\n        uploadData(file) {\n            if (!file) {\n                return;\n            }\n            let fr = new FileReader();\n            fr.addEventListener(\n                "load",\n                (e) => {\n                    let data = JSON.parse(e.target.result);\n                    this.fetchData(data);\n                }\n            );\n            fr.readAsText(file)\n        },\n        downloadData() {\n            if (!this.hits) {\n                return null;\n            }\n            const json = JSON.stringify(this.hits);\n            const blob = new Blob([json], { type: \'application/json\' });\n            const link = document.createElement(\'a\');\n            const date = new Date().toLocaleString(\'sv\').replace(\' \', \'_\').replaceAll(\'-\', \'_\').replaceAll(\':\', \'_\');\n            link.href = URL.createObjectURL(blob);\n            link.download = `${this.$STRINGS.APP_NAME}_${date}.json`;\n            link.click();\n            URL.revokeObjectURL(link.href);\n        },\n        resetProperties() {\n            this.ticket = "";\n            this.error = "";\n            this.mode = "";\n            this.hits = null;\n            this.selectedDatabases = 0;\n            this.tableMode = 0;\n        },\n        fetchData(data) {\n            this.resetProperties();\n            this.hits = [];\n            for (let result of data) {\n                this.hits.push(this.parseResults(result));\n            }\n        }\n    }\n};\n<\/script>\n\n<style scoped>\n::v-deep .v-app-bar-title__content {\n    text-overflow: revert !important;\n}\n</style>\n<style>\n.theme--light .panel-root .v-toolbar {\n    background-color: #454545 !important;\n}\n\n.theme--dark .panel-root .v-toolbar {\n    background-color: #1e1e1e !important;\n}\n</style>' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        226: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => o
            });
            var i = e(7537), a = e.n(i), r = e(3645), s = e.n(r)()(a());
            s.push([ t.id, "\n.structure-wrapper {\n    width: 400px;\n    height: 300px;\n    margin: 0 auto;\n}\n.theme--dark .structure-wrapper .v-tooltip__content {\n    background: rgba(97, 97, 97, 0.3);\n}\n/* @media only screen and (max-width: 600px) {\n    .structure-wrapper {\n        width: 300px;\n    }\n} */\n.structure-viewer {\n    width: 100%;\n    height: 100%;\n}\n.structure-viewer canvas {\n    border-radius: 2px;\n}\n.structure-panel {\n    position: relative;\n}\n.toolbar-panel {\n    display: inline-flex;\n    flex-direction: row;\n    position: absolute;\n    justify-content: center;\n    width: 100%;\n    bottom: 0;\n    z-index: 1;\n    left: 0;\n}\n.tmscore-panel {\n    position: absolute;\n    width: 100%;\n    top: 0;\n    left: 0;\n    z-index: 1;\n    font-family: monospace;\n    color: rgb(31, 119, 180);\n}\n.left-cell {\n    text-align: right;\n    width: 50%;\n}\n.right-cell {\n    text-align: left;\n    width: 50%;\n    padding-left: 0.3em;\n}\n", "", {
                version: 3,
                sources: [ "webpack://./frontend/StructureViewer.vue" ],
                names: [],
                mappings: ";AAujBA;IACA,YAAA;IACA,aAAA;IACA,cAAA;AACA;AAEA;IACA,iCAAA;AACA;AACA;;;;GAIA;AACA;IACA,WAAA;IACA,YAAA;AACA;AACA;IACA,kBAAA;AACA;AACA;IACA,kBAAA;AACA;AACA;IACA,oBAAA;IACA,mBAAA;IACA,kBAAA;IACA,uBAAA;IACA,WAAA;IACA,SAAA;IACA,UAAA;IACA,OAAA;AACA;AACA;IACA,kBAAA;IACA,WAAA;IACA,MAAA;IACA,OAAA;IACA,UAAA;IACA,sBAAA;IACA,wBAAA;AACA;AACA;IACA,iBAAA;IACA,UAAA;AACA;AACA;IACA,gBAAA;IACA,UAAA;IACA,mBAAA;AACA",
                sourcesContent: [ '<template>\n    <div class="structure-panel" v-if="\'tCa\' in alignment">\n        <div class="structure-wrapper" ref="structurepanel">\n            <v-tooltip open-delay="300" bottom attach=".structure-wrapper" background-color="transparent">\n                <template v-slot:activator="{ on }">\n                    <v-icon :light="isFullscreen" v-on="on" style="position: absolute; z-index: 999; right:0">{{ $MDI.HelpCircleOutline }}</v-icon>\n                </template>\n                <span>\n                    <dl style="text-align: center;">\n                        <dt>\n<svg xmlns="http://www.w3.org/2000/svg" xml:space="preserve" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linejoin:round;stroke-miterlimit:2" viewBox="0 0 32 32">\n<title>Left click</title>\n<path d="M25.6 5.8a5 5 0 0 0-5-4.8h-9.1a5 5 0 0 0-5.1 4.8v20.4a5 5 0 0 0 5 4.8h9.1a5 5 0 0 0 5.1-4.8V5.8Zm-1 9.5v10.9a4 4 0 0 1-4 3.8h-9.1a4 4 0 0 1-4-3.8V15.3h17ZM15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Zm1 0h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"/>\n<path id="left" d="M15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Z" style="fill:red"/>\n<path id="middle-inactive" d="M14.6 4h2.8v8h-2.8z"/>\n</svg>\n                        </dt>\n                        <dd>\n                            Rotate\n                        </dd>\n                        <dt>\n<svg xmlns="http://www.w3.org/2000/svg" xml:space="preserve" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linejoin:round;stroke-miterlimit:2" viewBox="0 0 32 32">\n<title>Right click</title>\n<path d="M25.6 5.8a5 5 0 0 0-5-4.8h-9.1a5 5 0 0 0-5.1 4.8v20.4a5 5 0 0 0 5 4.8h9.1a5 5 0 0 0 5.1-4.8V5.8Zm-1 9.5v10.9a4 4 0 0 1-4 3.8h-9.1a4 4 0 0 1-4-3.8V15.3h17ZM15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Zm1 0h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"/>\n<path id="right" d="M16.5 2h4a4 4 0 0 1 4 3.8v8.5h-8V2Z" style="fill:red"/>\n<path id="middle-inactive" d="M14.6 4h2.8v8h-2.8z"/>\n</svg>\n                        </dt>\n                        <dd>\n                            Pan\n                        </dd>\n                        <dt>\n<svg xmlns="http://www.w3.org/2000/svg" xml:space="preserve" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linejoin:round;stroke-miterlimit:2" viewBox="0 0 32 32">\n<title>Scroll wheel</title>\n<path d="M25.6 5.8a5 5 0 0 0-5-4.8h-9.1a5 5 0 0 0-5.1 4.8v20.4a5 5 0 0 0 5 4.8h9.1a5 5 0 0 0 5.1-4.8V5.8Zm-1 9.5v10.9a4 4 0 0 1-4 3.8h-9.1a4 4 0 0 1-4-3.8V15.3h17ZM15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Zm1 0h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"/>\n<path id="middle-active" d="M14.6 4h2.8v8h-2.8z" style="fill:red"/>\n</svg>\n                        </dt>\n                        <dd>\n                            Zoom\n                        </dd>\n                    </dl>\n                </span>\n            </v-tooltip>\n            <table v-if="tmAlignResults" class="tmscore-panel" v-bind="tmPanelBindings">\n                <tr>\n                    <td class="left-cell">TM-Score:</td>\n                    <td class="right-cell">{{ tmAlignResults.tmScore }}</td>\n                </tr>\n                <tr>\n                    <td class="left-cell">RMSD:</td>\n                    <td class="right-cell">{{ tmAlignResults.rmsd  }}</td>\n                </tr>\n            </table>\n            <div class="toolbar-panel">\n                <v-item-group class="v-btn-toggle" :light="isFullscreen">\n                <v-btn\n                    v-bind="tbButtonBindings"\n                    v-on:click="makePdb()"\n                    title="Save PDB"\n                >\n                    <v-icon v-bind="tbIconBindings">M19 3a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5c0-1.1.9-2 2-2h14Zm0 8v-.8c0-.7-.6-1.2-1.3-1.2h-2.4v6h2.4c.7 0 1.2-.5 1.2-1.2v-1c0-.4-.4-.8-.9-.8.5 0 1-.4 1-1Zm-9.7.5v-1c0-.8-.7-1.5-1.5-1.5H5.3v6h1.5v-2h1c.8 0 1.5-.7 1.5-1.5Zm5 2v-3c0-.8-.7-1.5-1.5-1.5h-2.5v6h2.5c.8 0 1.5-.7 1.5-1.5Zm3.4.3h-1.2v-1.2h1.2v1.2Zm-5.9-3.3v3h1v-3h-1Zm-5 0v1h1v-1h-1Zm11 .9h-1.3v-1.2h1.2v1.2Z</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Save PDB</span>\n                </v-btn>\n                <v-btn\n                    v-bind="tbButtonBindings"\n                    v-on:click="makeImage()"\n                    title="Save image"\n                >\n                    <v-icon v-bind="tbIconBindings">M19 3H5C3.9 3 3 3.9 3 5V19C3 20.1 3.9 21 5 21H19C20.1 21 21 20.1 21 19V5C21 3.9 20.1 3 19 3M9 11.5C9 12.3 8.3 13 7.5 13H6.5V15H5V9H7.5C8.3 9 9 9.7 9 10.5V11.5M14 15H12.5L11.5 12.5V15H10V9H11.5L12.5 11.5V9H14V15M19 10.5H16.5V13.5H17.5V12H19V13.7C19 14.4 18.5 15 17.7 15H16.4C15.6 15 15.1 14.3 15.1 13.7V10.4C15 9.7 15.5 9 16.3 9H17.6C18.4 9 18.9 9.7 18.9 10.3V10.5H19M6.5 10.5H7.5V11.5H6.5V10.5Z</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Save image</span>\n                </v-btn>\n                <v-btn\n                    v-bind="tbButtonBindings"\n                    v-on:click="cycleQueryView()"\n                    title="Toggle between the entire query structure and aligned region"\n                >\n                    <v-icon v-bind="tbIconBindings" style=\'color: #1E88E5;\' v-if="showQuery === 0">{{ ($LOCAL) ? $MDI.CircleHalf : "M12 12 V2 A10 10 0 0 0 3.858 17.806 Z" }}</v-icon>\n                    <v-icon v-bind="tbIconBindings" style=\'color: #1E88E5;\' v-else-if="!$LOCAL && showQuery === 1">M12 12 V2 A10 10 0 1 0 20.142 17.806 Z</v-icon>\n                    <v-icon v-bind="tbIconBindings" style=\'color: #1E88E5;\' v-else>{{ $MDI.Circle }}</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Toggle full query</span>\n              </v-btn>\n                <v-btn\n                    v-bind="tbButtonBindings"\n                    v-on:click="toggleFullTarget()"\n                    title="Toggle between the entire target structure and aligned region"\n                >\n                    <v-icon v-bind="tbIconBindings" style=\'color: #FFC107;\' v-if="showTarget == \'aligned\'">{{ $MDI.CircleHalf }}</v-icon>\n                    <v-icon v-bind="tbIconBindings" style=\'color: #FFC107;\' v-else>{{ $MDI.Circle }}</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Toggle full target</span>\n                </v-btn>\n                <v-btn\n                    v-bind="tbButtonBindings"\n                    v-on:click="toggleArrows()"\n                    title="Draw arrows between aligned residues"\n                >\n                    <v-icon v-bind="tbIconBindings" v-if="showArrows">{{ $MDI.ArrowRightCircle }}</v-icon>\n                    <v-icon v-bind="tbIconBindings" v-else>{{ $MDI.ArrowRightCircleOutline }}</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Toggle arrows</span>\n                </v-btn>\n                <v-btn\n                    v-bind="tbButtonBindings"\n                    v-on:click="resetView()"\n                    :input-value="\n                        selection != null\n                            && ((selection[0] != alignment.dbStartPos || selection[1] != alignment.dbEndPos)\n                            && (selection[0] != 1 || selection[1] != alignment.dbLen))"\n                    title="Reset the view to the original position and zoom level"\n                >\n                    <v-icon v-bind="tbIconBindings">{{ $MDI.Restore }}</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Reset view</span>\n                </v-btn>\n                <v-btn v-bind="tbButtonBindings"\n                    v-on:click="toggleFullscreen()"\n                    title="Enter fullscreen mode - press ESC to exit"\n                >\n                    <v-icon v-bind="tbIconBindings">{{ $MDI.Fullscreen }}</v-icon>\n                    <span v-if="isFullscreen">&nbsp;Fullscreen</span>\n                </v-btn>\n                </v-item-group>\n            </div>\n            <div class="structure-viewer" ref="viewport" />\n        </div>\n    </div>\n</template>\n\n<script>\nimport Panel from \'./Panel.vue\';\nimport { Shape, Stage, Selection, download, ColormakerRegistry, PdbWriter } from \'ngl\';\nimport { pulchra } from \'pulchra-wasm\';\nimport { tmalign, parse, parseMatrix } from \'tmalign-wasm\';\n\n\n// Create NGL arrows from array of ([X, Y, Z], [X, Y, Z]) pairs\nfunction createArrows(matches) {\n    const shape = new Shape(\'shape\')\n    for (let i = 0; i < matches.length; i++) {\n        const [a, b] = matches[i]\n        shape.addArrow(a, b, [0, 1, 1], 0.4)\n    }\n    return shape\n}\n\nconst oneToThree = {\n  "A":"ALA", "R":"ARG", "N":"ASN", "D":"ASP",\n  "C":"CYS", "E":"GLU", "Q":"GLN", "G":"GLY",\n  "H":"HIS", "I":"ILE", "L":"LEU", "K":"LYS",\n  "M":"MET", "F":"PHE", "P":"PRO", "S":"SER",\n  "T":"THR", "W":"TRP", "Y":"TYR", "V":"VAL",\n  "U":"SEC", "O":"PHL", "X":"XAA"\n};\n\n/**\n * Create a mock PDB from Ca data\n * Follows the spacing spec from https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM\n * Will have to change if/when swapping to fuller data\n */\nfunction mockPDB(ca, seq) {\n    const atoms = ca.split(\',\')\n    const pdb = new Array()\n    let j = 1\n    for (let i = 0; i < atoms.length; i += 3, j++) {\n        let [x, y, z] = atoms.slice(i, i + 3).map(element => parseFloat(element))\n        pdb.push(\n            \'ATOM  \'\n            + j.toString().padStart(5)\n            + \'  CA  \' + oneToThree[seq != "" && (atoms.length/3) == seq.length ? seq[i/3] : \'A\'] + \' A\'\n            + j.toString().padStart(4)\n            + \'    \'\n            + x.toString().padStart(8)\n            + y.toString().padStart(8)\n            + z.toString().padStart(8)\n            + \'  1.00  0.00           C  \'\n        )\n    }\n    return pdb.join(\'\\n\')\n}\n\n/* ------ The rotation matrix to rotate Chain_1 to Chain_2 ------ */\n/* m               t[m]        u[m][0]        u[m][1]        u[m][2] */\n/* 0     161.2708425765   0.0663961888  -0.6777150909  -0.7323208325 */\n/* 1     109.4205584665  -0.9559071424  -0.2536229340   0.1480437178 */\n/* 2      29.1924015422  -0.2860648199   0.6902011757  -0.6646722921 */\n/* Code for rotating Structure A from (x,y,z) to (X,Y,Z): */\n/* for(i=0; i<L; i++) */\n/* { */\n/*    X[i] = t[0] + u[0][0]*x[i] + u[0][1]*y[i] + u[0][2]*z[i]; */\n/*    Y[i] = t[1] + u[1][0]*x[i] + u[1][1]*y[i] + u[1][2]*z[i]; */\n/*    Z[i] = t[2] + u[2][0]*x[i] + u[2][1]*y[i] + u[2][2]*z[i]; */\n/* } */\nconst transformStructure = (structure, t, u) => {\n    structure.eachAtom(atom => {\n        const [x, y, z] = [atom.x, atom.y, atom.z]\n        atom.x = t[0] + u[0][0] * x + u[0][1] * y + u[0][2] * z\n        atom.y = t[1] + u[1][0] * x + u[1][1] * y + u[1][2] * z\n        atom.z = t[2] + u[2][0] * x + u[2][1] * y + u[2][2] * z\n    })\n    return structure\n}\n\n// Get XYZ coordinates of CA of a given residue\nconst xyz = (structure, resIndex) => {\n    var rp = structure.getResidueProxy()\n    var ap = structure.getAtomProxy()\n    rp.index = resIndex\n    ap.index = rp.getAtomIndexByName(\'CA\')\n    return [ap.x, ap.y, ap.z]\n}\n\n// Given an NGL AtomProxy, return the corresponding PDB line\nconst atomToPDBRow = (ap) => {\n    const { serial, atomname, resname, chainname, resno, inscode, x, y, z } = ap\n    return `ATOM  ${serial.toString().padStart(5)}${atomname.padStart(4)}  ${resname.padStart(3)} ${chainname.padStart(1)}${resno.toString().padStart(4)} ${inscode.padStart(1)}  ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}`\n}\n\n// Map 1-based indices in a selection to residue index/resno\nconst makeChainMap = (structure, sele) => {\n    let map = new Map()\n    structure.eachResidue(rp => {\n        map.set(rp.index + 1, { index: rp.index, resno: rp.resno });\n    }, new Selection(sele));\n    return map\n}\n\n// Generate a subsetted PDB file from a structure and selection\nconst makeSubPDB = (structure, sele) => {\n    let pdb = []\n    structure.eachAtom(ap => { pdb.push(atomToPDBRow(ap)) }, new Selection(sele))\n    return pdb.join(\'\\n\')\n}\n\nexport default {\n    components: { Panel },\n    data: () => ({\n        \'showTarget\': \'aligned\',\n        \'showQuery\': 0,\n        \'showArrows\': false,\n        \'selection\': null,\n        \'queryChain\': \'\',\n        \'qChainResMap\': null,\n        \'isFullscreen\': false,\n        \'tmAlignResults\': null,\n    }),\n    props: {\n        \'alignment\': Object,\n        \'queryFile\': String,\n        \'qColor\': { type: String, default: "white" },\n        \'tColor\': { type: String, default: "red" },\n        \'queryAlignedColor\': { type: String, default: "#1E88E5" },\n        \'queryUnalignedColor\': { type: String, default: "#A5CFF5" },\n        \'targetAlignedColor\': { type: String, default: "#FFC107" },\n        \'targetUnalignedColor\': { type: String, default: "#FFE699" },\n        \'qRepr\': { type: String, default: "cartoon" },\n        \'tRepr\': { type: String, default: "cartoon" },\n        \'bgColorLight\': { type: String, default: "white" },\n        \'bgColorDark\': { type: String, default: "#eee" },\n        \'queryMap\': { type: Array, default: null },\n        \'targetMap\': { type: Array, default: null },\n        \'hits\': { type: Object }\n    },\n    methods: {\n        // Parses two alignment strings, and saves matching residues\n        // Each match contains the index of the residue in the structure and a callback\n        // function to retrieve the residue\'s CA XYZ coordinates to allow retrieval\n        // before and after superposition (with updated coords)\n        saveMatchingResidues(aln1, aln2, str1, str2) {\n            if (aln1.length !== aln2.length) return\n            this.qMatches = []\n            this.tMatches = []\n            for (let i = 0; i < aln1.length; i++) {\n                if (aln1[i] === \'-\' || aln2[i] === \'-\') {\n                    continue;\n                }\n                // Make sure this residue actually exists in NGL structure representation\n                // e.g. d1b0ba starts with X, reported in alignment but removed by Pulchra\n                let qIdx = this.qChainResMap.get(this.queryMap[i]);\n                if (qIdx === undefined) {\n                    continue;\n                }\n                // Must be 0-based for xyz()\n                let tIdx = this.targetMap[i] - 1;\n                this.qMatches.push({ index: qIdx.index, xyz: () => xyz(str1, qIdx.index) })\n                this.tMatches.push({ index: tIdx, xyz: () => xyz(str2, tIdx) })\n            }\n        },\n        handleResize() {\n            if (!this.stage) return\n            this.stage.handleResize()\n        },\n        toggleFullscreen() {\n            if (!this.stage) return\n            this.stage.toggleFullscreen(this.$refs.structurepanel)\n        },\n        resetView() {\n            if (!this.stage) return\n            this.setSelection(this.showTarget)\n            this.stage.autoView(100)\n        },\n        toggleArrows() {\n            if (!this.stage || !this.arrowShape) return\n            this.showArrows = !this.showArrows\n        },\n        cycleQueryView() {\n            if (!this.stage)\n                return;\n            if (__LOCAL__) {\n                this.showQuery = (this.showQuery === 0) ? 1 : 0;\n            } else {\n                this.showQuery = (this.showQuery === 2) ? 0 : this.showQuery + 1;\n            }\n        },\n        toggleFullTarget() {\n            if (!this.stage) return\n            this.showTarget = this.showTarget === \'aligned\' ? \'full\' : \'aligned\'\n        },\n        setSelectionByRange(start, end) {\n            if (!this.targetRepr) return\n            this.targetRepr.setSelection(`${start}-${end}`)\n            this.stage.autoView(100)\n        },\n        setSelectionData(start, end) {\n            this.selection = [start, end]\n        },\n        setSelection(val) {\n            if (val === \'full\') this.setSelectionData(1, this.alignment.dbLen)\n            else this.setSelectionData(this.alignment.dbStartPos, this.alignment.dbEndPos)\n        },\n        setQuerySelection() {\n            this.queryRepr.setSelection(this.querySele)\n            this.stage.autoView(100)\n        },\n        // Update arrow shape on shape update\n        renderArrows() {\n            if (!this.stage) return\n            if (this.arrowShape) this.arrowShape.dispose()\n            let matches = new Array()\n            for (let i = 0; i < this.tMatches.length; i++) {\n                let qMatch = this.qMatches[i]\n                let tMatch = this.tMatches[i]\n                if (this.selection && !(tMatch.index >= this.selection[0] - 1 && tMatch.index < this.selection[1]))\n                    continue\n                matches.push([qMatch.xyz(), tMatch.xyz()])\n            }\n            this.arrowShape = this.stage.addComponentFromObject(createArrows(matches))\n            this.arrowShape.addRepresentation(\'buffer\')\n            this.arrowShape.setVisibility(this.showArrows)\n        },\n        makeImage() {\n            if (!this.stage) return\n            this.stage.viewer.setLight(undefined, undefined, undefined, 0.2)\n            this.stage.makeImage({\n                trim: true,\n                factor: (this.isFullscreen) ? 1 : 8,\n                antialias: true,\n                transparent: true,\n            }).then((blob) => {\n                this.stage.viewer.setLight(undefined, undefined, undefined, this.$vuetify.theme.dark ? 0.4 : 0.2)\n                download(blob, this.$route.params.ticket + \'-\' + this.alignment.target + ".png")\n            })\n        },\n        makePdb() {\n            if (!this.stage) return\n            if (!this.queryRepr) return\n            if (!this.targetRepr) return\n            let qPDB = new PdbWriter(this.queryRepr.repr.structure, { renumberSerial: false }).getData()\n            let tPDB = new PdbWriter(this.targetRepr.repr.structure, { renumberSerial: false }).getData()\n            qPDB = qPDB.split(\'\\n\').filter(line => line.startsWith(\'ATOM\')).join(\'\\n\')\n            tPDB = tPDB.split(\'\\n\').filter(line => line.startsWith(\'ATOM\')).join(\'\\n\')\n            let result =\n`TITLE     ${this.$route.params.ticket} - ${this.alignment.target}\nREMARK     This file was generated by the Foldseek webserver:\nREMARK       https://search.foldseek.com\nREMARK     Please cite:\nREMARK       https://doi.org/10.1101/2022.02.07.479398\nREMARK     Warning: Non C-alpha atoms might have been re-generated by PULCHRA,\nREMARK              if they are not present in the original PDB file.\nMODEL        1\n${qPDB}\nENDMDL\nMODEL        2\n${tPDB}\nENDMDL\nEND\n`\n            download(new Blob([result], { type: \'text/plain\' }), this.$route.params.ticket + \'-\' + this.alignment.target + ".pdb")\n        }\n    },\n    watch: {\n        \'showTarget\': function(val, _) {\n            this.setSelection(val)\n        },\n        \'showArrows\': function(val, _) {\n            if (!this.stage || !this.arrowShape) return\n            this.arrowShape.setVisibility(val)\n        },\n        \'selection\': function([start, end]) {\n            this.setSelectionByRange(start, end)\n            this.renderArrows()\n        },\n        \'showQuery\': function() {\n            if (!this.stage) return\n            this.setQuerySelection()\n        },\n        \'$route\': function() {}\n    },\n    computed: {\n        queryChainId: function() {\n            return this.queryChain.charCodeAt(0) - \'A\'.charCodeAt(0)\n        },\n        queryChainSele: function() {\n            return (this.queryChain) ? `(:${this.queryChain.toUpperCase()} OR :${this.queryChain.toLowerCase()})` : \'\';\n        },\n        querySubSele: function() {\n            if (!this.qChainResMap) {\n                return \'\';\n            }\n            let start = this.qChainResMap.get(this.alignment.qStartPos);\n            let end   = this.qChainResMap.get(this.alignment.qEndPos);\n            let sele  = `${start.resno}-${end.resno}`;\n            if (this.queryChain) {\n                sele = `${sele} AND ${this.queryChainSele}`;\n            }\n            return sele\n        },\n        querySele: function() {\n            if (this.showQuery == 0)\n                return this.querySubSele;\n            if (this.showQuery == 1)\n                return this.queryChainSele;\n            return \'\'\n        },\n        targetSele: function() {\n            if (!this.selection) return \'\'\n            return `${this.selection[0]}-${this.selection[1]}`;\n        },\n        tmPanelBindings: function() {\n            return (this.isFullscreen) ? { \'style\': \'margin-top: 10px; font-size: 2em; line-height: 2em\' } : {  }\n        },\n        tbIconBindings: function() {\n            return (this.isFullscreen) ? { \'right\': true } : {}\n        },\n        tbButtonBindings: function() {\n            return (this.isFullscreen) ? {\n                \'small\': false,\n                \'style\': \'margin-bottom: 15px;\',\n            } : {\n                \'small\': true,\n                \'style\': \'\'\n            }\n        }\n    },\n    beforeMount() {\n        let qChain = this.hits.query.header.match(/_([A-Z]+?)/m)\n        if (qChain) this.queryChain = qChain[0].replace(\'_\', \'\')\n    },\n    mounted() {\n        const bgColor = this.$vuetify.theme.dark ? this.bgColorDark : this.bgColorLight;\n        const ambientIntensity = this.$vuetify.theme.dark ? 0.4 : 0.2;\n        if (typeof(this.alignment.tCa) == "undefined")\n            return;\n        this.stage = new Stage(this.$refs.viewport,\n            {\n                backgroundColor: bgColor,\n                ambientIntensity: ambientIntensity,\n                clipNear: -1000,\n                clipFar: 1000,\n                fogFar: 1000,\n                fogNear: -1000,\n                quality: \'high\'\n            })\n        \n        Promise.all([\n            (this.$LOCAL)\n                ? pulchra(mockPDB(this.hits.query.qCa, this.hits.query.sequence))\n                : this.$axios.get("api/result/" + this.$route.params.ticket + \'/query\'),\n            pulchra(mockPDB(this.alignment.tCa, this.alignment.tSeq))\n        ]).then(([qResponse, tPdb]) => {\n            // Sanitize PDB in case of lines with too few characters\n            qResponse = (__LOCAL__) ? qResponse : qResponse.data;\n            let data = \'\';\n            for (let line of qResponse.split(\'\\n\')) {\n                let numCols = Math.max(0, 80 - line.length);\n                let newLine = line + \' \'.repeat(numCols) + \'\\n\';\n                data += newLine\n            }\n            qResponse = data;\n            return [qResponse, tPdb];\n        }).then(([qResponse, tPdb]) => {\n            Promise.all([\n                this.stage.loadFile(new Blob([qResponse], { type: \'text/plain\' }), {ext: \'pdb\', firstModelOnly: true}),\n                this.stage.loadFile(new Blob([tPdb], { type: \'text/plain\' }), {ext: \'pdb\', firstModelOnly: true}),\n            ])\n            .then(([query, target]) => {\n                // fix query if its course-grained (i.e. non-complete backbone)\n                // foldseek doesn\'t care about the oxygen, but NGL does\n                if (query.structure.getAtomProxy().isCg()) {\n                    return new Promise((resolve, reject) => {\n                        pulchra(qResponse)\n                        .then((queryPdb) => {\n                            this.stage.loadFile(new Blob([queryPdb], { type: \'text/plain\' }), {ext: \'pdb\', firstModelOnly: true})\n                            .then((query) => { resolve([query, target]) })\n                            .catch(reject)\n                        }).catch(reject)\n                    })\n                } else {\n                    return [query, target];\n                }\n            })\n            .then(([query, target]) => {\n                // Map 1-based indices to residue index/resno; only need for query structure\n                // Use queryChainSele to make all selections based on actual query chain\n                this.qChainResMap = makeChainMap(query.structure, this.queryChainSele)\n                this.saveMatchingResidues(this.alignment.qAln, this.alignment.dbAln, query.structure, target.structure)\n\n                // Generate colorschemes for query/target based on alignment\n                this.querySchemeId = ColormakerRegistry.addSelectionScheme([\n                    [this.queryAlignedColor, this.querySubSele],\n                    [this.queryUnalignedColor, "*"],\n                ], "_queryScheme")\n                this.targetSchemeId = ColormakerRegistry.addSelectionScheme([\n                    [this.targetAlignedColor, `${this.alignment.dbStartPos}-${this.alignment.dbEndPos}`],\n                    [this.targetUnalignedColor, "*"]\n                ], "_targetScheme")\n\n                // Generate subsetted PDBs for TM-align\n                let qSubPdb = makeSubPDB(query.structure, this.querySubSele)\n                let tSubPdb = makeSubPDB(target.structure, `${this.alignment.dbStartPos}-${this.alignment.dbEndPos}`)\n                let alnFasta = `>target\\n${this.alignment.dbAln}\\n\\n>query\\n${this.alignment.qAln}`\n\n                // Re-align target to query using TM-align for better superposition\n                // Target 1st since TM-align generates superposition matrix for 1st structure\n                tmalign(tSubPdb, qSubPdb, alnFasta).then(out => {\n                    this.tmAlignResults = parse(out.output)\n                    let { t, u } = parseMatrix(out.matrix)\n                    transformStructure(target.structure, t, u)\n                    this.queryRepr = query.addRepresentation(this.qRepr, {color: this.querySchemeId})\n                    this.targetRepr = target.addRepresentation(this.tRepr, {color: this.targetSchemeId})\n                }).then(() => {\n                    this.setSelection(this.showTarget)\n                    this.setQuerySelection()\n                    this.stage.autoView()\n                })\n            })\n        })\n        window.addEventListener(\'resize\', this.handleResize)\n        this.stage.signals.fullscreenChanged.add((isFullscreen) => {\n            if (isFullscreen) {\n                this.stage.viewer.setBackground(\'#ffffff\')\n                this.stage.viewer.setLight(undefined, undefined, undefined, 0.2)\n                this.isFullscreen = true\n            } else {\n                this.stage.viewer.setBackground(bgColor)\n                this.stage.viewer.setLight(undefined, undefined, undefined, ambientIntensity)\n                this.isFullscreen = false\n            }\n        })\n    },\n    beforeDestroy() {\n        if (typeof(this.stage) == \'undefined\')\n            return\n        this.stage.dispose() \n        window.removeEventListener(\'resize\', this.handleResize)\n    }\n}\n<\/script>\n\n<style>\n.structure-wrapper {\n    width: 400px;\n    height: 300px;\n    margin: 0 auto;\n}\n\n.theme--dark .structure-wrapper .v-tooltip__content {\n    background: rgba(97, 97, 97, 0.3);\n}\n/* @media only screen and (max-width: 600px) {\n    .structure-wrapper {\n        width: 300px;\n    }\n} */\n.structure-viewer {\n    width: 100%;\n    height: 100%;\n}\n.structure-viewer canvas {\n    border-radius: 2px;\n}\n.structure-panel {\n    position: relative;\n}\n.toolbar-panel {\n    display: inline-flex;\n    flex-direction: row;\n    position: absolute;\n    justify-content: center;\n    width: 100%;\n    bottom: 0;\n    z-index: 1;\n    left: 0;\n}\n.tmscore-panel {\n    position: absolute;\n    width: 100%;\n    top: 0;\n    left: 0;\n    z-index: 1;\n    font-family: monospace;\n    color: rgb(31, 119, 180);\n}\n.left-cell {\n    text-align: right;\n    width: 50%;\n}\n.right-cell {\n    text-align: left;\n    width: 50%;\n    padding-left: 0.3em;\n}\n</style>\n' ],
                sourceRoot: ""
            } ]);
            const o = s;
        },
        9010: (t, n, e) => {
            var i = e(7537), a = e(3645), r = e(1667), s = e(7204), o = e(1464), l = a(i), A = r(s), c = r(o);
            l.push([ t.id, "@font-face{font-family:InconsolataClustal;src:url(" + A + "),url(" + c + ')}.hide{display:none}.db{border-left:5px solid #000}@media print,screen and (max-width: 599px){small.ticket{display:inline-block;line-height:.9}}.result-table a.anchor{display:block;position:relative;top:-125px;visibility:hidden}.result-table a:not([href]){color:#333}.result-table a:not([href]):not([href]):hover{text-decoration:none}.result-table td,.result-table th{padding:0 6px;text-align:left}.result-table .hit.active{background:#f9f9f9}.result-table tbody:hover td[rowspan],.result-table tbody tr:hover{background:#eee}.result-table .alignment-action{text-align:center;word-wrap:normal}.theme--dark .result-table a:not([href]){color:#eee}.theme--dark .result-table .hit.active{background:#333}.theme--dark .result-table tbody:hover td[rowspan],.theme--dark .result-table tbody tr:hover{background:#333}@media print,screen and (min-width: 961px){.result-table{table-layout:fixed;border-collapse:collapse;width:100%}.result-table th.wide-1{width:15%}.result-table th.wide-2{width:30%}.result-table th.wide-3{width:45%}.result-table th.thin{width:6.5% !important;white-space:nowrap}.result-table td.thin{white-space:nowrap}.result-table .long{overflow:hidden;word-break:keep-all;text-overflow:ellipsis;white-space:nowrap}}@media print{.result-table .alignment-action{display:none}}@media screen and (max-width: 960px){.result-table{width:100%}.result-table .long{height:100% !important;white-space:normal !important;min-height:48px}.result-table .hits{min-width:300px}.result-table tbody td a{min-width:100px}.result-table tbody td.graphical div.ruler{margin:10px 0}.result-table thead{display:none}.result-table tfoot th{border:0;display:inherit}.result-table tr{box-shadow:0 2px 3px rgba(0,0,0,.1),0 0 0 1px rgba(0,0,0,.1);max-width:100%;position:relative;display:block;padding:.5em}.result-table tr td{border:0;display:inherit}.result-table tr td:last-child{border-bottom:0}.result-table tr:not(:last-child){margin-bottom:1rem}.result-table tr:not(.is-selected){background:inherit}.result-table tr:not(.is-selected):hover{background-color:inherit}.result-table tr.detail{margin-top:-1rem}.result-table tr:not(.detail):not(.is-empty):not(.table-footer) td{display:flex;border-bottom:1px solid #eee;flex-direction:row}.result-table tr:not(.detail):not(.is-empty):not(.table-footer) td:last-child{border-bottom:0}.result-table tr:not(.detail):not(.is-empty):not(.table-footer) td:before{content:attr(data-label);font-weight:600;margin-right:auto;padding-right:.5em;word-break:keep-all;flex:1;white-space:nowrap}.result-table tbody td a,.result-table tbody td span{flex:2;margin-left:auto;text-align:right;word-wrap:anywhere}}.alignment{position:absolute;left:4px;right:4px;z-index:999}.alignment .residues{font-family:InconsolataClustal,Inconsolata,Consolas,Menlo,Monaco,"Cascadia Mono","Segoe UI Mono","Roboto Mono","Oxygen Mono","Ubuntu Monospace","Source Code Pro","Fira Mono","Droid Sans Mono","Courier New",monospace;white-space:pre}.theme--dark .alignment .residues{color:#fff}.clear-button{font:14px sans-serif;cursor:pointer}', "", {
                version: 3,
                sources: [ "webpack://./frontend/ResultView.vue" ],
                names: [],
                mappings: "AAoRA,WACA,8BAAA,CACA,mFAAA,CAIA,MACI,YAAA,CAGJ,IACI,0BAAA,CAGJ,2CACA,aACI,oBAAA,CACA,cAAA,CAAA,CAKA,uBACI,aAAA,CACA,iBAAA,CACA,UAAA,CACA,iBAAA,CAGJ,4BACI,UAAA,CACA,8CACI,oBAAA,CAIR,kCACI,aAAA,CACA,eAAA,CAGJ,0BACI,kBAAA,CAGJ,mEACI,eAAA,CAGJ,gCACI,iBAAA,CACA,gBAAA,CAOA,yCACI,UAAA,CAGJ,uCACI,eAAA,CAGJ,6FACI,eAAA,CAKZ,2CACI,cACI,kBAAA,CACA,wBAAA,CACA,UAAA,CACA,wBACI,SAAA,CAEJ,wBACI,SAAA,CAGJ,wBACI,SAAA,CAEJ,sBACI,qBAAA,CACA,kBAAA,CAEJ,sBACI,kBAAA,CAEJ,oBACI,eAAA,CACA,mBAAA,CACA,sBAAA,CACA,kBAAA,CAAA,CAKZ,aACI,gCACI,YAAA,CAAA,CAIR,qCACI,cACI,UAAA,CACA,oBACI,sBAAA,CACA,6BAAA,CACA,eAAA,CAGJ,oBACI,eAAA,CAGJ,yBACI,eAAA,CAGJ,2CACI,aAAA,CAGJ,oBACI,YAAA,CAGJ,uBACI,QAAA,CACA,eAAA,CAGJ,iBACI,4DAAA,CACA,cAAA,CACA,iBAAA,CACA,aAAA,CACA,YAAA,CAGJ,oBACI,QAAA,CACA,eAAA,CAGJ,+BACI,eAAA,CAEJ,kCACI,kBAAA,CAEJ,mCACI,kBAAA,CAEJ,yCACI,wBAAA,CAEJ,wBACI,gBAAA,CAGJ,mEACI,YAAA,CACA,4BAAA,CACA,kBAAA,CAEA,8EACI,eAAA,CAGR,0EACI,wBAAA,CACA,eAAA,CACA,iBAAA,CACA,kBAAA,CACA,mBAAA,CACA,MAAA,CACA,kBAAA,CAGJ,qDACI,MAAA,CACA,gBAAA,CACA,gBAAA,CACA,kBAAA,CAAA,CAKZ,WACI,iBAAA,CACA,QAAA,CACA,SAAA,CACA,WAAA,CAEA,qBACI,uNAAA,CACA,eAAA,CAIA,kCACI,UAAA,CAKZ,cACI,oBAAA,CACA,cAAA",
                sourcesContent: [ '@import "_variables.scss";\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n@font-face {\nfont-family: InconsolataClustal;\nsrc: url(assets/InconsolataClustal2.woff2),\n     url(assets/InconsolataClustal2.woff);\n}\n\n.hide {\n    display: none;\n}\n\n.db {\n    border-left: 5px solid black;\n}\n\n@media print, screen and (max-width: 599px) {\nsmall.ticket {\n    display: inline-block;\n    line-height: 0.9;\n}\n}\n\n.result-table {\n    a.anchor {\n        display: block;\n        position: relative;\n        top: -125px;\n        visibility: hidden;\n    }\n\n    a:not([href]) {\n        color: #333;\n        &:not([href]):hover {\n            text-decoration: none;\n        }\n    }\n\n    td, th {\n        padding: 0 6px;\n        text-align: left;\n    }\n\n    .hit.active {\n        background: #f9f9f9;\n    }\n\n    tbody:hover td[rowspan], tbody tr:hover {\n        background: #eee;\n    }\n\n    .alignment-action {\n        text-align: center;\n        word-wrap: normal;\n    }\n}\n\n\n.theme--dark {\n    .result-table {\n        a:not([href])  {\n            color: #eee;\n        }\n\n        .hit.active {\n            background: #333;\n        }\n\n        tbody:hover td[rowspan], tbody tr:hover {\n            background: #333;\n        }\n    }\n}\n\n@media print, screen and (min-width: 961px) {\n    .result-table {\n        table-layout: fixed;\n        border-collapse: collapse;\n        width: 100%;\n        th.wide-1 {\n            width: 15%;\n        }\n        th.wide-2 {\n            width: 30%;\n        }\n\n        th.wide-3 {\n            width: 45%;\n        }\n        th.thin {\n            width: 6.5% !important;\n            white-space: nowrap;\n        }\n        td.thin {\n            white-space: nowrap;\n        }\n        .long {\n            overflow: hidden;\n            word-break: keep-all;\n            text-overflow: ellipsis;\n            white-space: nowrap;\n        }\n    }\n}\n\n@media print {\n    .result-table .alignment-action {\n        display: none;\n    }\n}\n\n@media screen and (max-width: 960px) {\n    .result-table {\n        width: 100%;\n        .long {\n            height: 100% !important;\n            white-space: normal !important;\n            min-height: 48px;\n        }\n\n        .hits {\n            min-width: 300px;\n        }\n\n        tbody td a {\n            min-width: 100px;\n        }\n\n        tbody td.graphical div.ruler {\n            margin: 10px 0;\n        }\n\n        thead {\n            display: none;\n        }\n\n        tfoot th {\n            border: 0;\n            display: inherit;\n        }\n\n        tr {\n            box-shadow: 0 2px 3px rgba(0, 0, 0, 0.1), 0 0 0 1px rgba(0, 0, 0, 0.1);\n            max-width: 100%;\n            position: relative;\n            display: block;\n            padding: 0.5em;\n        }\n\n        tr td {\n            border: 0;\n            display: inherit;\n        }\n\n        tr td:last-child {\n            border-bottom: 0;\n        }\n        tr:not(:last-child) {\n            margin-bottom: 1rem;\n        }\n        tr:not(.is-selected) {\n            background: inherit;\n        }\n        tr:not(.is-selected):hover {\n            background-color: inherit;\n        }\n        tr.detail {\n            margin-top: -1rem;\n        }\n\n        tr:not(.detail):not(.is-empty):not(.table-footer) td {\n            display: flex;\n            border-bottom: 1px solid #eee;\n            flex-direction: row;\n\n            &:last-child {\n                border-bottom: 0;\n            }\n        }\n        tr:not(.detail):not(.is-empty):not(.table-footer) td:before {\n            content: attr(data-label);\n            font-weight: 600;\n            margin-right: auto;\n            padding-right: 0.5em;\n            word-break: keep-all;\n            flex: 1;\n            white-space: nowrap;\n        }\n\n        tbody td a, tbody td span {\n            flex: 2;\n            margin-left: auto;\n            text-align: right;\n            word-wrap: anywhere;\n        }\n    }\n}\n\n.alignment {\n    position:absolute;\n    left:4px;\n    right:4px;\n    z-index: 999;\n\n    .residues {\n        font-family: InconsolataClustal, Inconsolata, Consolas, Menlo, Monaco, "Cascadia Mono", "Segoe UI Mono", "Roboto Mono", "Oxygen Mono", "Ubuntu Monospace", "Source Code Pro", "Fira Mono", "Droid Sans Mono", "Courier New", monospace;\n        white-space: pre;\n    }\n\n    .theme--dark & {\n        .residues {\n            color: #fff;\n        }\n    }\n}\n\n.clear-button {\n    font: 14px sans-serif;\n    cursor: pointer;\n}\n\n\n' ],
                sourceRoot: ""
            } ]), t.exports = l;
        },
        5385: (t, n, e) => {
            var i = e(7537), a = e(3645)(i);
            a.push([ t.id, ".ruler[data-v-2b7861b2]{position:relative;width:100%;height:10px;border-top:1px solid #333}.tick-label[data-v-2b7861b2]{position:absolute;word-wrap:normal;font-size:9px;word-break:keep-all;line-height:1em;margin-top:7px;width:50px;margin-left:-25px;text-align:center;font-weight:bold}.tick-label-top[data-v-2b7861b2]{margin-top:-15px}.query[data-v-2b7861b2]{position:absolute;top:0;bottom:0;margin-top:-5px;--chevron-width: 5px;height:10px}.chevron-start[data-v-2b7861b2]{position:absolute;left:0;bottom:0;top:0;width:5px;clip-path:polygon(0 0, var(--chevron-width) 0, var(--chevron-width) 100%, 0 100%, var(--chevron-width) 50%)}.query.reversed .chevron-start[data-v-2b7861b2]{clip-path:polygon(var(--chevron-width) 0, 0 50%, var(--chevron-width) 100%)}.chevron-mid[data-v-2b7861b2]{position:absolute;left:5px;right:5px;bottom:0;top:0}.chevron-end[data-v-2b7861b2]{position:absolute;right:0;bottom:0;top:0;width:5px;clip-path:polygon(0 0, var(--chevron-width) 50%, 0 100%)}.query.reversed .chevron-end[data-v-2b7861b2]{clip-path:polygon(0 0, var(--chevron-width) 0, 0 50%, var(--chevron-width) 100%, 0 100%);clip-path:polygon()}.theme--dark .ruler[data-v-2b7861b2]{border-color:#aaa}", "", {
                version: 3,
                sources: [ "webpack://./frontend/Ruler.vue" ],
                names: [],
                mappings: "AAwDA,wBACE,iBAAA,CACA,UAAA,CACA,WAAA,CACA,yBAAA,CAGF,6BACE,iBAAA,CACA,gBAAA,CACA,aAAA,CACA,mBAAA,CACA,eAAA,CACA,cAAA,CACA,UAAA,CACA,iBAAA,CACA,iBAAA,CACA,gBAAA,CAGF,iCACE,gBAAA,CAGF,wBACE,iBAAA,CACA,KAAA,CACA,QAAA,CACA,eAAA,CACA,oBAAA,CACA,WAAA,CAGF,gCACE,iBAAA,CACA,MAAA,CACA,QAAA,CACA,KAAA,CACA,SAAA,CACA,2GAAA,CAGF,gDACE,2EAAA,CAGF,8BACE,iBAAA,CACA,QAAA,CACA,SAAA,CACA,QAAA,CACA,KAAA,CAGF,8BACE,iBAAA,CACA,OAAA,CACA,QAAA,CACA,KAAA,CACA,SAAA,CACA,wDAAA,CAEF,8CACE,wFAAA,CACA,mBAAA,CAIE,qCACE,iBAAA",
                sourcesContent: [ '@import "_variables.scss";\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n.ruler {\n  position: relative;\n  width: 100%;\n  height: 10px;\n  border-top: 1px solid #333;\n}\n\n.tick-label {\n  position: absolute;\n  word-wrap: normal;\n  font-size: 9px;\n  word-break: keep-all;\n  line-height: 1em;\n  margin-top: 7px;\n  width: 50px;\n  margin-left: -25px;\n  text-align: center;\n  font-weight: bold;\n}\n\n.tick-label-top {\n  margin-top: -15px;\n}\n\n.query {\n  position: absolute;\n  top: 0;\n  bottom: 0;\n  margin-top: -5px;\n  --chevron-width: 5px;\n  height: 10px;\n}\n\n.chevron-start {\n  position: absolute;\n  left:0;\n  bottom:0;\n  top:0;\n  width:5px;\n  clip-path: polygon(0 0, var(--chevron-width) 0, var(--chevron-width) 100%, 0 100%, var(--chevron-width) 50%);\n}\n\n.query.reversed .chevron-start {\n  clip-path: polygon(var(--chevron-width) 0, 0 50%, var(--chevron-width) 100%);\n}\n\n.chevron-mid {\n  position: absolute;\n  left:5px;\n  right:5px;\n  bottom:0;\n  top:0;\n}\n\n.chevron-end {\n  position: absolute;\n  right:0;\n  bottom:0;\n  top:0;\n  width:5px;\n  clip-path: polygon(0 0, var(--chevron-width) 50%, 0 100%);\n}\n.query.reversed .chevron-end {\n  clip-path: polygon(0 0, var(--chevron-width) 0, 0 50%, var(--chevron-width) 100%, 0 100%);\n  clip-path: polygon()\n}\n\n.theme--dark {\n    .ruler {\n      border-color: #aaa;\n    }\n}\n' ],
                sourceRoot: ""
            } ]), t.exports = a;
        },
        654: (t, n, e) => {
            var i = e(9837);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("4fa110d4", i, !1, {});
        },
        603: (t, n, e) => {
            var i = e(5426);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("59383ee7", i, !1, {});
        },
        2530: (t, n, e) => {
            var i = e(6696);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("4a805097", i, !1, {});
        },
        4449: (t, n, e) => {
            var i = e(8260);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("28e700a6", i, !1, {});
        },
        9146: (t, n, e) => {
            var i = e(4569);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("5d44b975", i, !1, {});
        },
        2556: (t, n, e) => {
            var i = e(864);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("0a2d9f56", i, !1, {});
        },
        8973: (t, n, e) => {
            var i = e(8742);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("77ba9bdc", i, !1, {});
        },
        6608: (t, n, e) => {
            var i = e(226);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("1147822a", i, !1, {});
        },
        5264: (t, n, e) => {
            var i = e(9010);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("122feea2", i, !1, {});
        },
        5941: (t, n, e) => {
            var i = e(5385);
            i.__esModule && (i = i.default), "string" == typeof i && (i = [ [ t.id, i, "" ] ]), 
            i.locals && (t.exports = i.locals);
            (0, e(5346).Z)("6d831950", i, !1, {});
        },
        917: (t, n, e) => {
            "use strict";
            e.d(n, {
                Z: () => m
            });
            var i = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return e("div", {
                    class: [ "panel-root", null != t.elevation ? "elevation-" + t.elevation : null ]
                }, [ t.$slots.header || t.header ? e("v-toolbar", {
                    attrs: {
                        text: "",
                        dense: "",
                        dark: ""
                    }
                }, [ t.collapsible ? e("v-btn", {
                    staticStyle: {
                        "margin-top": "0",
                        "margin-left": "-15px"
                    },
                    attrs: {
                        icon: "",
                        plain: "",
                        "aria-expanded": t.isCollapsed ? "false" : "true",
                        "aria-controls": t.uuid
                    },
                    on: {
                        click: function(n) {
                            t.isCollapsed = !t.isCollapsed;
                        }
                    }
                }, [ t.isCollapsed ? e("v-icon", [ t._v("\n                " + t._s(t.$MDI.PlusBox) + "\n            ") ]) : e("v-icon", [ t._v("\n                " + t._s(t.$MDI.MinusBox) + "\n            ") ]) ], 1) : t._e(), t._v(" "), e("span", {
                    staticClass: "text-h6 align-end"
                }, [ t.$slots.header ? t._t("header") : [ t._v(t._s(t.header)) ] ], 2), t._v(" "), e("v-spacer"), t._v(" "), t._t("toolbar-extra") ], 2) : t._e(), t._v(" "), t.isCollapsed ? t._e() : e("v-card", {
                    class: [ "panel", {
                        "d-flex": t.flex
                    }, {
                        "force-fill-height": t.fillHeight
                    } ],
                    attrs: {
                        rounded: "0",
                        id: t.uuid
                    }
                }, [ t.$slots.desc ? e("v-card-text", {
                    staticClass: "subheading justify"
                }, [ t._t("desc") ], 2) : t._e(), t._v(" "), t.$slots.content ? e("v-card-text", {
                    class: [ "panel-content", "justify", {
                        "d-flex": t.flex
                    } ]
                }, [ t._t("content") ], 2) : t._e() ], 1) ], 1);
            };
            i._withStripped = !0;
            var a = 0;
            const r = {
                name: "panel",
                props: {
                    header: {
                        default: "",
                        type: String
                    },
                    fillHeight: {
                        default: !1,
                        type: Boolean
                    },
                    collapsible: {
                        default: !1,
                        type: Boolean
                    },
                    collapsed: {
                        default: !1,
                        type: Boolean
                    },
                    flex: {
                        default: !0,
                        type: Boolean
                    },
                    elevation: {
                        default: null,
                        type: Number
                    }
                },
                data: function() {
                    return {
                        isCollapsed: this.collapsed
                    };
                },
                beforeCreate: function() {
                    this.uuid = "panel-" + a.toString(), a += 1;
                }
            };
            e(9146);
            var s = e(1900), o = e(3453), l = e.n(o), A = e(5934), c = e(5893), d = e(5255), u = e(4786), h = e(2515), p = e(9567), g = (0, 
            s.Z)(r, i, [], !1, null, "0d9b5935", null);
            l()(g, {
                VBtn: A.Z,
                VCard: c.Z,
                VCardText: d.ZB,
                VIcon: u.Z,
                VSpacer: h.Z,
                VToolbar: p.Z
            }), g.options.__file = "frontend/Panel.vue";
            const m = g.exports;
        },
        8992: (t, n, e) => {
            "use strict";
            e.r(n), e.d(n, {
                default: () => x
            });
            var i = function() {
                var t = this, n = t.$createElement, e = t._self._c || n;
                return "tCa" in t.alignment ? e("div", {
                    staticClass: "structure-panel"
                }, [ e("div", {
                    ref: "structurepanel",
                    staticClass: "structure-wrapper"
                }, [ e("v-tooltip", {
                    attrs: {
                        "open-delay": "300",
                        bottom: "",
                        attach: ".structure-wrapper",
                        "background-color": "transparent"
                    },
                    scopedSlots: t._u([ {
                        key: "activator",
                        fn: function(n) {
                            var i = n.on;
                            return [ e("v-icon", t._g({
                                staticStyle: {
                                    position: "absolute",
                                    "z-index": "999",
                                    right: "0"
                                },
                                attrs: {
                                    light: t.isFullscreen
                                }
                            }, i), [ t._v(t._s(t.$MDI.HelpCircleOutline)) ]) ];
                        }
                    } ], null, !1, 1827187420)
                }, [ t._v(" "), e("span", [ e("dl", {
                    staticStyle: {
                        "text-align": "center"
                    }
                }, [ e("dt", [ e("svg", {
                    staticStyle: {
                        "fill-rule": "evenodd",
                        "clip-rule": "evenodd",
                        "stroke-linejoin": "round",
                        "stroke-miterlimit": "2"
                    },
                    attrs: {
                        xmlns: "http://www.w3.org/2000/svg",
                        "xml:space": "preserve",
                        viewBox: "0 0 32 32"
                    }
                }, [ e("title", [ t._v("Left click") ]), t._v(" "), e("path", {
                    attrs: {
                        d: "M25.6 5.8a5 5 0 0 0-5-4.8h-9.1a5 5 0 0 0-5.1 4.8v20.4a5 5 0 0 0 5 4.8h9.1a5 5 0 0 0 5.1-4.8V5.8Zm-1 9.5v10.9a4 4 0 0 1-4 3.8h-9.1a4 4 0 0 1-4-3.8V15.3h17ZM15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Zm1 0h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"
                    }
                }), t._v(" "), e("path", {
                    staticStyle: {
                        fill: "red"
                    },
                    attrs: {
                        id: "left",
                        d: "M15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Z"
                    }
                }), t._v(" "), e("path", {
                    attrs: {
                        id: "middle-inactive",
                        d: "M14.6 4h2.8v8h-2.8z"
                    }
                }) ]) ]), t._v(" "), e("dd", [ t._v("\n                            Rotate\n                        ") ]), t._v(" "), e("dt", [ e("svg", {
                    staticStyle: {
                        "fill-rule": "evenodd",
                        "clip-rule": "evenodd",
                        "stroke-linejoin": "round",
                        "stroke-miterlimit": "2"
                    },
                    attrs: {
                        xmlns: "http://www.w3.org/2000/svg",
                        "xml:space": "preserve",
                        viewBox: "0 0 32 32"
                    }
                }, [ e("title", [ t._v("Right click") ]), t._v(" "), e("path", {
                    attrs: {
                        d: "M25.6 5.8a5 5 0 0 0-5-4.8h-9.1a5 5 0 0 0-5.1 4.8v20.4a5 5 0 0 0 5 4.8h9.1a5 5 0 0 0 5.1-4.8V5.8Zm-1 9.5v10.9a4 4 0 0 1-4 3.8h-9.1a4 4 0 0 1-4-3.8V15.3h17ZM15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Zm1 0h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"
                    }
                }), t._v(" "), e("path", {
                    staticStyle: {
                        fill: "red"
                    },
                    attrs: {
                        id: "right",
                        d: "M16.5 2h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"
                    }
                }), t._v(" "), e("path", {
                    attrs: {
                        id: "middle-inactive",
                        d: "M14.6 4h2.8v8h-2.8z"
                    }
                }) ]) ]), t._v(" "), e("dd", [ t._v("\n                            Pan\n                        ") ]), t._v(" "), e("dt", [ e("svg", {
                    staticStyle: {
                        "fill-rule": "evenodd",
                        "clip-rule": "evenodd",
                        "stroke-linejoin": "round",
                        "stroke-miterlimit": "2"
                    },
                    attrs: {
                        xmlns: "http://www.w3.org/2000/svg",
                        "xml:space": "preserve",
                        viewBox: "0 0 32 32"
                    }
                }, [ e("title", [ t._v("Scroll wheel") ]), t._v(" "), e("path", {
                    attrs: {
                        d: "M25.6 5.8a5 5 0 0 0-5-4.8h-9.1a5 5 0 0 0-5.1 4.8v20.4a5 5 0 0 0 5 4.8h9.1a5 5 0 0 0 5.1-4.8V5.8Zm-1 9.5v10.9a4 4 0 0 1-4 3.8h-9.1a4 4 0 0 1-4-3.8V15.3h17ZM15.5 2v12.3h-8V5.8a4 4 0 0 1 4-3.8h4Zm1 0h4a4 4 0 0 1 4 3.8v8.5h-8V2Z"
                    }
                }), t._v(" "), e("path", {
                    staticStyle: {
                        fill: "red"
                    },
                    attrs: {
                        id: "middle-active",
                        d: "M14.6 4h2.8v8h-2.8z"
                    }
                }) ]) ]), t._v(" "), e("dd", [ t._v("\n                            Zoom\n                        ") ]) ]) ]) ]), t._v(" "), t.tmAlignResults ? e("table", t._b({
                    staticClass: "tmscore-panel"
                }, "table", t.tmPanelBindings, !1), [ e("tr", [ e("td", {
                    staticClass: "left-cell"
                }, [ t._v("TM-Score:") ]), t._v(" "), e("td", {
                    staticClass: "right-cell"
                }, [ t._v(t._s(t.tmAlignResults.tmScore)) ]) ]), t._v(" "), e("tr", [ e("td", {
                    staticClass: "left-cell"
                }, [ t._v("RMSD:") ]), t._v(" "), e("td", {
                    staticClass: "right-cell"
                }, [ t._v(t._s(t.tmAlignResults.rmsd)) ]) ]) ]) : t._e(), t._v(" "), e("div", {
                    staticClass: "toolbar-panel"
                }, [ e("v-item-group", {
                    staticClass: "v-btn-toggle",
                    attrs: {
                        light: t.isFullscreen
                    }
                }, [ e("v-btn", t._b({
                    attrs: {
                        title: "Save PDB"
                    },
                    on: {
                        click: function(n) {
                            return t.makePdb();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ e("v-icon", t._b({}, "v-icon", t.tbIconBindings, !1), [ t._v("M19 3a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5c0-1.1.9-2 2-2h14Zm0 8v-.8c0-.7-.6-1.2-1.3-1.2h-2.4v6h2.4c.7 0 1.2-.5 1.2-1.2v-1c0-.4-.4-.8-.9-.8.5 0 1-.4 1-1Zm-9.7.5v-1c0-.8-.7-1.5-1.5-1.5H5.3v6h1.5v-2h1c.8 0 1.5-.7 1.5-1.5Zm5 2v-3c0-.8-.7-1.5-1.5-1.5h-2.5v6h2.5c.8 0 1.5-.7 1.5-1.5Zm3.4.3h-1.2v-1.2h1.2v1.2Zm-5.9-3.3v3h1v-3h-1Zm-5 0v1h1v-1h-1Zm11 .9h-1.3v-1.2h1.2v1.2Z") ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Save PDB") ]) : t._e() ], 1), t._v(" "), e("v-btn", t._b({
                    attrs: {
                        title: "Save image"
                    },
                    on: {
                        click: function(n) {
                            return t.makeImage();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ e("v-icon", t._b({}, "v-icon", t.tbIconBindings, !1), [ t._v("M19 3H5C3.9 3 3 3.9 3 5V19C3 20.1 3.9 21 5 21H19C20.1 21 21 20.1 21 19V5C21 3.9 20.1 3 19 3M9 11.5C9 12.3 8.3 13 7.5 13H6.5V15H5V9H7.5C8.3 9 9 9.7 9 10.5V11.5M14 15H12.5L11.5 12.5V15H10V9H11.5L12.5 11.5V9H14V15M19 10.5H16.5V13.5H17.5V12H19V13.7C19 14.4 18.5 15 17.7 15H16.4C15.6 15 15.1 14.3 15.1 13.7V10.4C15 9.7 15.5 9 16.3 9H17.6C18.4 9 18.9 9.7 18.9 10.3V10.5H19M6.5 10.5H7.5V11.5H6.5V10.5Z") ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Save image") ]) : t._e() ], 1), t._v(" "), e("v-btn", t._b({
                    attrs: {
                        title: "Toggle between the entire query structure and aligned region"
                    },
                    on: {
                        click: function(n) {
                            return t.cycleQueryView();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ 0 === t.showQuery ? e("v-icon", t._b({
                    staticStyle: {
                        color: "#1E88E5"
                    }
                }, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$LOCAL ? t.$MDI.CircleHalf : "M12 12 V2 A10 10 0 0 0 3.858 17.806 Z")) ]) : t.$LOCAL || 1 !== t.showQuery ? e("v-icon", t._b({
                    staticStyle: {
                        color: "#1E88E5"
                    }
                }, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.Circle)) ]) : e("v-icon", t._b({
                    staticStyle: {
                        color: "#1E88E5"
                    }
                }, "v-icon", t.tbIconBindings, !1), [ t._v("M12 12 V2 A10 10 0 1 0 20.142 17.806 Z") ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Toggle full query") ]) : t._e() ], 1), t._v(" "), e("v-btn", t._b({
                    attrs: {
                        title: "Toggle between the entire target structure and aligned region"
                    },
                    on: {
                        click: function(n) {
                            return t.toggleFullTarget();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ "aligned" == t.showTarget ? e("v-icon", t._b({
                    staticStyle: {
                        color: "#FFC107"
                    }
                }, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.CircleHalf)) ]) : e("v-icon", t._b({
                    staticStyle: {
                        color: "#FFC107"
                    }
                }, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.Circle)) ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Toggle full target") ]) : t._e() ], 1), t._v(" "), e("v-btn", t._b({
                    attrs: {
                        title: "Draw arrows between aligned residues"
                    },
                    on: {
                        click: function(n) {
                            return t.toggleArrows();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ t.showArrows ? e("v-icon", t._b({}, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.ArrowRightCircle)) ]) : e("v-icon", t._b({}, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.ArrowRightCircleOutline)) ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Toggle arrows") ]) : t._e() ], 1), t._v(" "), e("v-btn", t._b({
                    attrs: {
                        "input-value": !(null == t.selection || t.selection[0] == t.alignment.dbStartPos && t.selection[1] == t.alignment.dbEndPos || 1 == t.selection[0] && t.selection[1] == t.alignment.dbLen),
                        title: "Reset the view to the original position and zoom level"
                    },
                    on: {
                        click: function(n) {
                            return t.resetView();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ e("v-icon", t._b({}, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.Restore)) ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Reset view") ]) : t._e() ], 1), t._v(" "), e("v-btn", t._b({
                    attrs: {
                        title: "Enter fullscreen mode - press ESC to exit"
                    },
                    on: {
                        click: function(n) {
                            return t.toggleFullscreen();
                        }
                    }
                }, "v-btn", t.tbButtonBindings, !1), [ e("v-icon", t._b({}, "v-icon", t.tbIconBindings, !1), [ t._v(t._s(t.$MDI.Fullscreen)) ]), t._v(" "), t.isFullscreen ? e("span", [ t._v(" Fullscreen") ]) : t._e() ], 1) ], 1) ], 1), t._v(" "), e("div", {
                    ref: "viewport",
                    staticClass: "structure-viewer"
                }) ], 1) ]) : t._e();
            };
            i._withStripped = !0;
            var a = e(8152), r = e(917), s = e(8197), o = e(7895), l = e(1434);
            function A(t, n) {
                var e = "undefined" != typeof Symbol && t[Symbol.iterator] || t["@@iterator"];
                if (!e) {
                    if (Array.isArray(t) || (e = function(t, n) {
                        if (!t) return;
                        if ("string" == typeof t) return c(t, n);
                        var e = Object.prototype.toString.call(t).slice(8, -1);
                        "Object" === e && t.constructor && (e = t.constructor.name);
                        if ("Map" === e || "Set" === e) return Array.from(t);
                        if ("Arguments" === e || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(e)) return c(t, n);
                    }(t)) || n && t && "number" == typeof t.length) {
                        e && (t = e);
                        var i = 0, a = function() {};
                        return {
                            s: a,
                            n: function() {
                                return i >= t.length ? {
                                    done: !0
                                } : {
                                    done: !1,
                                    value: t[i++]
                                };
                            },
                            e: function(t) {
                                throw t;
                            },
                            f: a
                        };
                    }
                    throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.");
                }
                var r, s = !0, o = !1;
                return {
                    s: function() {
                        e = e.call(t);
                    },
                    n: function() {
                        var t = e.next();
                        return s = t.done, t;
                    },
                    e: function(t) {
                        o = !0, r = t;
                    },
                    f: function() {
                        try {
                            s || null == e.return || e.return();
                        } finally {
                            if (o) throw r;
                        }
                    }
                };
            }
            function c(t, n) {
                (null == n || n > t.length) && (n = t.length);
                for (var e = 0, i = new Array(n); e < n; e++) i[e] = t[e];
                return i;
            }
            var d = {
                A: "ALA",
                R: "ARG",
                N: "ASN",
                D: "ASP",
                C: "CYS",
                E: "GLU",
                Q: "GLN",
                G: "GLY",
                H: "HIS",
                I: "ILE",
                L: "LEU",
                K: "LYS",
                M: "MET",
                F: "PHE",
                P: "PRO",
                S: "SER",
                T: "THR",
                W: "TRP",
                Y: "TYR",
                V: "VAL",
                U: "SEC",
                O: "PHL",
                X: "XAA"
            };
            function u(t, n) {
                for (var e = t.split(","), i = new Array, r = 1, s = 0; s < e.length; s += 3, r++) {
                    var o = e.slice(s, s + 3).map((function(t) {
                        return parseFloat(t);
                    })), l = (0, a.Z)(o, 3), A = l[0], c = l[1], u = l[2];
                    i.push("ATOM  " + r.toString().padStart(5) + "  CA  " + d["" != n && e.length / 3 == n.length ? n[s / 3] : "A"] + " A" + r.toString().padStart(4) + "    " + A.toString().padStart(8) + c.toString().padStart(8) + u.toString().padStart(8) + "  1.00  0.00           C  ");
                }
                return i.join("\n");
            }
            var h = function(t, n) {
                var e = t.getResidueProxy(), i = t.getAtomProxy();
                return e.index = n, i.index = e.getAtomIndexByName("CA"), [ i.x, i.y, i.z ];
            }, p = function(t, n) {
                var e = [];
                return t.eachAtom((function(t) {
                    e.push(function(t) {
                        var n = t.serial, e = t.atomname, i = t.resname, a = t.chainname, r = t.resno, s = t.inscode, o = t.x, l = t.y, A = t.z;
                        return "ATOM  ".concat(n.toString().padStart(5)).concat(e.padStart(4), "  ").concat(i.padStart(3), " ").concat(a.padStart(1)).concat(r.toString().padStart(4), " ").concat(s.padStart(1), "  ").concat(o.toFixed(3).padStart(8)).concat(l.toFixed(3).padStart(8)).concat(A.toFixed(3).padStart(8));
                    }(t));
                }), new s.Y1(n)), e.join("\n");
            };
            const g = {
                components: {
                    Panel: r.Z
                },
                data: function() {
                    return {
                        showTarget: "aligned",
                        showQuery: 0,
                        showArrows: !1,
                        selection: null,
                        queryChain: "",
                        qChainResMap: null,
                        isFullscreen: !1,
                        tmAlignResults: null
                    };
                },
                props: {
                    alignment: Object,
                    queryFile: String,
                    qColor: {
                        type: String,
                        default: "white"
                    },
                    tColor: {
                        type: String,
                        default: "red"
                    },
                    queryAlignedColor: {
                        type: String,
                        default: "#1E88E5"
                    },
                    queryUnalignedColor: {
                        type: String,
                        default: "#A5CFF5"
                    },
                    targetAlignedColor: {
                        type: String,
                        default: "#FFC107"
                    },
                    targetUnalignedColor: {
                        type: String,
                        default: "#FFE699"
                    },
                    qRepr: {
                        type: String,
                        default: "cartoon"
                    },
                    tRepr: {
                        type: String,
                        default: "cartoon"
                    },
                    bgColorLight: {
                        type: String,
                        default: "white"
                    },
                    bgColorDark: {
                        type: String,
                        default: "#eee"
                    },
                    queryMap: {
                        type: Array,
                        default: null
                    },
                    targetMap: {
                        type: Array,
                        default: null
                    },
                    hits: {
                        type: Object
                    }
                },
                methods: {
                    saveMatchingResidues: function(t, n, e, i) {
                        var a = this;
                        if (t.length === n.length) {
                            this.qMatches = [], this.tMatches = [];
                            for (var r = function() {
                                if ("-" === t[s] || "-" === n[s]) return 0;
                                var r = a.qChainResMap.get(a.queryMap[s]);
                                if (void 0 === r) return 0;
                                var o = a.targetMap[s] - 1;
                                a.qMatches.push({
                                    index: r.index,
                                    xyz: function() {
                                        return h(e, r.index);
                                    }
                                }), a.tMatches.push({
                                    index: o,
                                    xyz: function() {
                                        return h(i, o);
                                    }
                                });
                            }, s = 0; s < t.length; s++) r();
                        }
                    },
                    handleResize: function() {
                        this.stage && this.stage.handleResize();
                    },
                    toggleFullscreen: function() {
                        this.stage && this.stage.toggleFullscreen(this.$refs.structurepanel);
                    },
                    resetView: function() {
                        this.stage && (this.setSelection(this.showTarget), this.stage.autoView(100));
                    },
                    toggleArrows: function() {
                        this.stage && this.arrowShape && (this.showArrows = !this.showArrows);
                    },
                    cycleQueryView: function() {
                        this.stage && (this.showQuery = 0 === this.showQuery ? 1 : 0);
                    },
                    toggleFullTarget: function() {
                        this.stage && (this.showTarget = "aligned" === this.showTarget ? "full" : "aligned");
                    },
                    setSelectionByRange: function(t, n) {
                        this.targetRepr && (this.targetRepr.setSelection("".concat(t, "-").concat(n)), this.stage.autoView(100));
                    },
                    setSelectionData: function(t, n) {
                        this.selection = [ t, n ];
                    },
                    setSelection: function(t) {
                        "full" === t ? this.setSelectionData(1, this.alignment.dbLen) : this.setSelectionData(this.alignment.dbStartPos, this.alignment.dbEndPos);
                    },
                    setQuerySelection: function() {
                        this.queryRepr.setSelection(this.querySele), this.stage.autoView(100);
                    },
                    renderArrows: function() {
                        if (this.stage) {
                            this.arrowShape && this.arrowShape.dispose();
                            for (var t = new Array, n = 0; n < this.tMatches.length; n++) {
                                var e = this.qMatches[n], i = this.tMatches[n];
                                (!this.selection || i.index >= this.selection[0] - 1 && i.index < this.selection[1]) && t.push([ e.xyz(), i.xyz() ]);
                            }
                            this.arrowShape = this.stage.addComponentFromObject(function(t) {
                                for (var n = new s.bn("shape"), e = 0; e < t.length; e++) {
                                    var i = (0, a.Z)(t[e], 2), r = i[0], o = i[1];
                                    n.addArrow(r, o, [ 0, 1, 1 ], .4);
                                }
                                return n;
                            }(t)), this.arrowShape.addRepresentation("buffer"), this.arrowShape.setVisibility(this.showArrows);
                        }
                    },
                    makeImage: function() {
                        var t = this;
                        this.stage && (this.stage.viewer.setLight(void 0, void 0, void 0, .2), this.stage.makeImage({
                            trim: !0,
                            factor: this.isFullscreen ? 1 : 8,
                            antialias: !0,
                            transparent: !0
                        }).then((function(n) {
                            t.stage.viewer.setLight(void 0, void 0, void 0, t.$vuetify.theme.dark ? .4 : .2), 
                            (0, s.LR)(n, t.$route.params.ticket + "-" + t.alignment.target + ".png");
                        })));
                    },
                    makePdb: function() {
                        if (this.stage && this.queryRepr && this.targetRepr) {
                            var t = new s.p8(this.queryRepr.repr.structure, {
                                renumberSerial: !1
                            }).getData(), n = new s.p8(this.targetRepr.repr.structure, {
                                renumberSerial: !1
                            }).getData();
                            t = t.split("\n").filter((function(t) {
                                return t.startsWith("ATOM");
                            })).join("\n"), n = n.split("\n").filter((function(t) {
                                return t.startsWith("ATOM");
                            })).join("\n");
                            var e = "TITLE     ".concat(this.$route.params.ticket, " - ").concat(this.alignment.target, "\nREMARK     This file was generated by the Foldseek webserver:\nREMARK       https://search.foldseek.com\nREMARK     Please cite:\nREMARK       https://doi.org/10.1101/2022.02.07.479398\nREMARK     Warning: Non C-alpha atoms might have been re-generated by PULCHRA,\nREMARK              if they are not present in the original PDB file.\nMODEL        1\n").concat(t, "\nENDMDL\nMODEL        2\n").concat(n, "\nENDMDL\nEND\n");
                            (0, s.LR)(new Blob([ e ], {
                                type: "text/plain"
                            }), this.$route.params.ticket + "-" + this.alignment.target + ".pdb");
                        }
                    }
                },
                watch: {
                    showTarget: function(t, n) {
                        this.setSelection(t);
                    },
                    showArrows: function(t, n) {
                        this.stage && this.arrowShape && this.arrowShape.setVisibility(t);
                    },
                    selection: function(t) {
                        var n = (0, a.Z)(t, 2), e = n[0], i = n[1];
                        this.setSelectionByRange(e, i), this.renderArrows();
                    },
                    showQuery: function() {
                        this.stage && this.setQuerySelection();
                    },
                    $route: function() {}
                },
                computed: {
                    queryChainId: function() {
                        return this.queryChain.charCodeAt(0) - "A".charCodeAt(0);
                    },
                    queryChainSele: function() {
                        return this.queryChain ? "(:".concat(this.queryChain.toUpperCase(), " OR :").concat(this.queryChain.toLowerCase(), ")") : "";
                    },
                    querySubSele: function() {
                        if (!this.qChainResMap) return "";
                        var t = this.qChainResMap.get(this.alignment.qStartPos), n = this.qChainResMap.get(this.alignment.qEndPos), e = "".concat(t.resno, "-").concat(n.resno);
                        return this.queryChain && (e = "".concat(e, " AND ").concat(this.queryChainSele)), 
                        e;
                    },
                    querySele: function() {
                        return 0 == this.showQuery ? this.querySubSele : 1 == this.showQuery ? this.queryChainSele : "";
                    },
                    targetSele: function() {
                        return this.selection ? "".concat(this.selection[0], "-").concat(this.selection[1]) : "";
                    },
                    tmPanelBindings: function() {
                        return this.isFullscreen ? {
                            style: "margin-top: 10px; font-size: 2em; line-height: 2em"
                        } : {};
                    },
                    tbIconBindings: function() {
                        return this.isFullscreen ? {
                            right: !0
                        } : {};
                    },
                    tbButtonBindings: function() {
                        return this.isFullscreen ? {
                            small: !1,
                            style: "margin-bottom: 15px;"
                        } : {
                            small: !0,
                            style: ""
                        };
                    }
                },
                beforeMount: function() {
                    var t = this.hits.query.header.match(/_([A-Z]+?)/m);
                    t && (this.queryChain = t[0].replace("_", ""));
                },
                mounted: function() {
                    var t = this, n = this.$vuetify.theme.dark ? this.bgColorDark : this.bgColorLight, e = this.$vuetify.theme.dark ? .4 : .2;
                    void 0 !== this.alignment.tCa && (this.stage = new s.Hf(this.$refs.viewport, {
                        backgroundColor: n,
                        ambientIntensity: e,
                        clipNear: -1e3,
                        clipFar: 1e3,
                        fogFar: 1e3,
                        fogNear: -1e3,
                        quality: "high"
                    }), Promise.all([ this.$LOCAL ? (0, o.n)(u(this.hits.query.qCa, this.hits.query.sequence)) : this.$axios.get("api/result/" + this.$route.params.ticket + "/query"), (0, 
                    o.n)(u(this.alignment.tCa, this.alignment.tSeq)) ]).then((function(t) {
                        var n, e = (0, a.Z)(t, 2), i = e[0], r = e[1], s = "", o = A(i.split("\n"));
                        try {
                            for (o.s(); !(n = o.n()).done; ) {
                                var l = n.value, c = Math.max(0, 80 - l.length);
                                s += l + " ".repeat(c) + "\n";
                            }
                        } catch (t) {
                            o.e(t);
                        } finally {
                            o.f();
                        }
                        return [ i = s, r ];
                    })).then((function(n) {
                        var e = (0, a.Z)(n, 2), i = e[0], r = e[1];
                        Promise.all([ t.stage.loadFile(new Blob([ i ], {
                            type: "text/plain"
                        }), {
                            ext: "pdb",
                            firstModelOnly: !0
                        }), t.stage.loadFile(new Blob([ r ], {
                            type: "text/plain"
                        }), {
                            ext: "pdb",
                            firstModelOnly: !0
                        }) ]).then((function(n) {
                            var e = (0, a.Z)(n, 2), r = e[0], s = e[1];
                            return r.structure.getAtomProxy().isCg() ? new Promise((function(n, e) {
                                (0, o.n)(i).then((function(i) {
                                    t.stage.loadFile(new Blob([ i ], {
                                        type: "text/plain"
                                    }), {
                                        ext: "pdb",
                                        firstModelOnly: !0
                                    }).then((function(t) {
                                        n([ t, s ]);
                                    })).catch(e);
                                })).catch(e);
                            })) : [ r, s ];
                        })).then((function(n) {
                            var e, i, r, o = (0, a.Z)(n, 2), A = o[0], c = o[1];
                            t.qChainResMap = (e = A.structure, i = t.queryChainSele, r = new Map, e.eachResidue((function(t) {
                                r.set(t.index + 1, {
                                    index: t.index,
                                    resno: t.resno
                                });
                            }), new s.Y1(i)), r), t.saveMatchingResidues(t.alignment.qAln, t.alignment.dbAln, A.structure, c.structure), 
                            t.querySchemeId = s.Ub.addSelectionScheme([ [ t.queryAlignedColor, t.querySubSele ], [ t.queryUnalignedColor, "*" ] ], "_queryScheme"), 
                            t.targetSchemeId = s.Ub.addSelectionScheme([ [ t.targetAlignedColor, "".concat(t.alignment.dbStartPos, "-").concat(t.alignment.dbEndPos) ], [ t.targetUnalignedColor, "*" ] ], "_targetScheme");
                            var d = p(A.structure, t.querySubSele), u = p(c.structure, "".concat(t.alignment.dbStartPos, "-").concat(t.alignment.dbEndPos)), h = ">target\n".concat(t.alignment.dbAln, "\n\n>query\n").concat(t.alignment.qAln);
                            (0, l.Mb)(u, d, h).then((function(n) {
                                t.tmAlignResults = (0, l.Qc)(n.output);
                                var e = (0, l.im)(n.matrix), i = e.t, a = e.u;
                                !function(t, n, e) {
                                    t.eachAtom((function(t) {
                                        var i = [ t.x, t.y, t.z ], a = i[0], r = i[1], s = i[2];
                                        t.x = n[0] + e[0][0] * a + e[0][1] * r + e[0][2] * s, t.y = n[1] + e[1][0] * a + e[1][1] * r + e[1][2] * s, 
                                        t.z = n[2] + e[2][0] * a + e[2][1] * r + e[2][2] * s;
                                    }));
                                }(c.structure, i, a), t.queryRepr = A.addRepresentation(t.qRepr, {
                                    color: t.querySchemeId
                                }), t.targetRepr = c.addRepresentation(t.tRepr, {
                                    color: t.targetSchemeId
                                });
                            })).then((function() {
                                t.setSelection(t.showTarget), t.setQuerySelection(), t.stage.autoView();
                            }));
                        }));
                    })), window.addEventListener("resize", this.handleResize), this.stage.signals.fullscreenChanged.add((function(i) {
                        i ? (t.stage.viewer.setBackground("#ffffff"), t.stage.viewer.setLight(void 0, void 0, void 0, .2), 
                        t.isFullscreen = !0) : (t.stage.viewer.setBackground(n), t.stage.viewer.setLight(void 0, void 0, void 0, e), 
                        t.isFullscreen = !1);
                    })));
                },
                beforeDestroy: function() {
                    void 0 !== this.stage && (this.stage.dispose(), window.removeEventListener("resize", this.handleResize));
                }
            };
            e(6608);
            var m = e(1900), v = e(3453), f = e.n(v), b = e(5934), C = e(4786), M = e(7309), y = e(1562), w = (0, 
            m.Z)(g, i, [], !1, null, null, null);
            f()(w, {
                VBtn: b.Z,
                VIcon: C.Z,
                VItemGroup: M.Z,
                VTooltip: y.Z
            }), w.options.__file = "frontend/StructureViewer.vue";
            const x = w.exports;
        },
        1464: t => {
            "use strict";
            t.exports = "data:font/woff;base64,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";
        },
        7204: t => {
            "use strict";
            t.exports = "data:font/woff2;base64,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";
        },
        4833: (t, n, e) => {
            "use strict";
            t.exports = e.p + "142d6904f2305dd1cce7.png";
        },
        5904: (t, n, e) => {
            "use strict";
            t.exports = e.p + "9f772eefe8d08175ff5d.png";
        },
        6617: (t, n, e) => {
            "use strict";
            t.exports = e.p + "70a4e8d38900d34cea12.png";
        },
        7969: (t, n, e) => {
            "use strict";
            t.exports = e.p + "53d2a61fad6a2df4af57.png";
        },
        5515: (t, n, e) => {
            "use strict";
            t.exports = e.p + "fb5cfc3806f721f541ad.png";
        },
        4484: (t, n, e) => {
            "use strict";
            t.exports = e.p + "cb013a3d1b5f9a2c78e2.png";
        },
        7940: (t, n, e) => {
            "use strict";
            t.exports = e.p + "753a136eb8e7d5534788.png";
        },
        7018: (t, n, e) => {
            "use strict";
            t.exports = e.p + "9c4014f243b1c404a691.png";
        },
        2881: (t, n, e) => {
            "use strict";
            t.exports = e.p + "ea7744ed67559f380a81.svg";
        },
        42: (t, n, e) => {
            "use strict";
            t.exports = e.p + "a1ee785acc7f8c1bf4ac.png";
        },
        901: (t, n, e) => {
            "use strict";
            t.exports = e.p + "bfc0aaa54b3fd8130101.png";
        }
    }, e = {};
    function i(t) {
        var a = e[t];
        if (void 0 !== a) return a.exports;
        var r = e[t] = {
            id: t,
            loaded: !1,
            exports: {}
        };
        return n[t].call(r.exports, r, r.exports, i), r.loaded = !0, r.exports;
    }
    i.m = n, t = [], i.O = (n, e, a, r) => {
        if (!e) {
            var s = 1 / 0;
            for (c = 0; c < t.length; c++) {
                for (var [e, a, r] = t[c], o = !0, l = 0; l < e.length; l++) (!1 & r || s >= r) && Object.keys(i.O).every((t => i.O[t](e[l]))) ? e.splice(l--, 1) : (o = !1, 
                r < s && (s = r));
                if (o) {
                    t.splice(c--, 1);
                    var A = a();
                    void 0 !== A && (n = A);
                }
            }
            return n;
        }
        r = r || 0;
        for (var c = t.length; c > 0 && t[c - 1][2] > r; c--) t[c] = t[c - 1];
        t[c] = [ e, a, r ];
    }, i.n = t => {
        var n = t && t.__esModule ? () => t.default : () => t;
        return i.d(n, {
            a: n
        }), n;
    }, i.d = (t, n) => {
        for (var e in n) i.o(n, e) && !i.o(t, e) && Object.defineProperty(t, e, {
            enumerable: !0,
            get: n[e]
        });
    }, i.g = function() {
        if ("object" == typeof globalThis) return globalThis;
        try {
            return this || new Function("return this")();
        } catch (t) {
            if ("object" == typeof window) return window;
        }
    }(), i.o = (t, n) => Object.prototype.hasOwnProperty.call(t, n), i.r = t => {
        "undefined" != typeof Symbol && Symbol.toStringTag && Object.defineProperty(t, Symbol.toStringTag, {
            value: "Module"
        }), Object.defineProperty(t, "__esModule", {
            value: !0
        });
    }, i.nmd = t => (t.paths = [], t.children || (t.children = []), t), i.p = "/", (() => {
        i.b = document.baseURI || self.location.href;
        var t = {
            179: 0
        };
        i.O.j = n => 0 === t[n];
        var n = (n, e) => {
            var a, r, [s, o, l] = e, A = 0;
            if (s.some((n => 0 !== t[n]))) {
                for (a in o) i.o(o, a) && (i.m[a] = o[a]);
                if (l) var c = l(i);
            }
            for (n && n(e); A < s.length; A++) r = s[A], i.o(t, r) && t[r] && t[r][0](), t[r] = 0;
            return i.O(c);
        }, e = self.webpackChunkmmseqs_app = self.webpackChunkmmseqs_app || [];
        e.forEach(n.bind(null, 0)), e.push = n.bind(null, e.push.bind(e));
    })();
    var a = i.O(void 0, [ 736 ], (() => i(8217)));
    a = i.O(a);
})();
//# sourceMappingURL=main.js.map