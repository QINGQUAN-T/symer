# -*- coding: utf-8 -*-
"""W7 图表: 论文图 3/4/5/6/7/8/9/10 统一渲染（v4 稿图号）。
数据源: results/*.json（图级源数据, 由 w3_figs_data / w4_lone_spreader /
w6_network_ledger / w6_failure_boundary / w7_measure / w7_sim_archives /
w7_network_ci 生成）。面板标号 (a)(b)(c) 统一右上角。
输出: results/figs/fig{3,4,5,6,7,8,9,10}-*.pdf"""
import json
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300,
    "font.family": "serif", "font.size": 9,
    "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.4,
    "mathtext.fontset": "cm",
})
OI = {"blue": "#0072B2", "verm": "#D55E00", "green": "#009E73",
      "purple": "#CC79A7", "orange": "#E69E00", "black": "#000000"}
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGD = f'{ROOT}/results/figs'

os.makedirs(FIGD, exist_ok=True)
sim = json.load(open(f'{ROOT}/results/w7_simarchive_results.json'))
net = json.load(open(f'{ROOT}/results/w7_network_ci.json'))

def _load(name):
    try:
        return json.load(open(f'{ROOT}/results/{name}'))
    except FileNotFoundError:
        print(f'[w7_figs] 缺少 results/{name} —— 先运行对应数据脚本; 跳过该图')
        return None

def tag(ax, s):
    ax.set_title(s, loc='left', fontsize=10)

# ---------------- Fig 10 (v4: simulated archives) ----------------
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.0))
TH = 2.0
ax = axes[0]
for mu, col, mk in [(0.5, OI['blue'], 'o'), (1.37, OI['verm'], 's'), (2.0, OI['green'], '^')]:
    xs, ys, lo, hi = [], [], [], []
    for c in sim['cells']:
        if c['mu'] == mu and c['theta'] == TH and c['alpha'] > 0:
            xs.append(c['alpha']); ys.append(c['ahat_med'])
            lo.append(c['ahat_med'] - c['ahat_lo']); hi.append(c['ahat_hi'] - c['ahat_med'])
    ax.errorbar(xs, ys, yerr=[lo, hi], marker=mk, ms=4, lw=1, color=col,
                capsize=2, label=fr'$\mu={mu:g}$')
ax.plot([2e-3, 1.2e-1], [2e-3, 1.2e-1], '--', lw=0.8, color='0.6')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('true $\\alpha$')
ax.set_ylabel(r'inferred $\hat\alpha$ (median, IQR)')
ax.set_ylim(6e-4, 2.5e-1)
ax.legend(loc='upper left', frameon=False, fontsize=7)
ax.text(0.97, 0.06, 'no false positives at $\\alpha=0$ (all 9 cells)',
        transform=ax.transAxes, fontsize=6.8, color='0.35', va='bottom', ha='right')
tag(ax, '(a)')

ax = axes[1]
for mu, col, mk in [(0.5, OI['blue'], 'o'), (1.37, OI['verm'], 's'), (2.0, OI['green'], '^')]:
    xs, ys = [], []
    for c in sim['cells']:
        if c['mu'] == mu and c['theta'] == TH and c['alpha'] > 0 and c['relbias_med'] is not None:
            xs.append(c['alpha']); ys.append(c['relbias_med'] * 100)
    ax.plot(xs, ys, marker=mk, ms=4, lw=1, color=col, label=fr'$\mu={mu:g}$')
ax.axhline(0, color='0.6', ls='--', lw=0.8)
ax.set_xscale('log')
ax.set_xlabel('true $\\alpha$')
ax.set_ylabel('relative bias of $\\hat\\alpha$ (%)  at $\\theta=2$')
ax.set_ylim(-88, 78)
ax.legend(loc='center right', frameon=False, fontsize=7)
tag(ax, '(b)')

ax = axes[2]
deg = sim['degradation']
labs, vals, err = [], [], []
for c in deg:
    if c['drop'] > 0:
        lab = f"drop {int(c['drop']*100)}%"
    else:
        lab = {1.0: 'baseline', 0.9: 'span 90%', 0.75: 'span 75%'}[c['censor']]
    labs.append(lab)
    vals.append(c['ahat_med']); err.append(max(c['ahat_hi'] - c['ahat_med'],
                                               c['ahat_med'] - c['ahat_lo']))
cols = [OI['blue']] + ['#56B4E9'] * (len(vals) - 1)
ax.bar(range(len(vals)), vals, yerr=err, capsize=2.5, color=cols, width=0.62)
ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, rotation=30, ha='right',
                                                    fontsize=7.5)
ax.set_ylim(0, 0.031)
ax.set_ylabel(r'median $\hat\alpha$')
tag(ax, '(c)')
fig.tight_layout()
fig.savefig(f'{FIGD}/fig10-simarchive.pdf', bbox_inches='tight')
print('fig10-simarchive saved')

# ---------------- Fig 7 (v4: network CI) ----------------
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
ax = axes[0]
pts = net['exp']
zz = np.array([p['k'] for p in pts], float)
dv = np.array([p['dev_mean'] for p in pts], float)
er = np.array([p['dev_sd'] for p in pts], float) * 1.96
ax.errorbar(zz, dv, yerr=er, fmt='o-', ms=4, lw=1, color=OI['blue'],
            capsize=2, label='simulation (10 graphs $\\times$ 3 seeds)')
fit = net['fit']
sl, ic = np.polyfit(np.log(zz), np.log(dv), 1)
xs = np.geomspace(zz.min() * 0.85, zz.max() * 1.15, 50)
ax.plot(xs, np.exp(ic) * xs ** sl, '--', lw=0.9, color=OI['orange'],
        label=fr'fit $\langle k\rangle^{{-{fit["exponent"]:.2f}}}$ '
              fr'[{-fit["hi"]:.2f}, {-fit["lo"]:.2f}]')
ax.set_xlabel(r'mean degree $\langle k\rangle$')
ax.set_ylabel(r'ledger deviation (%)')
ax.set_ylim(0, 16)
ax.legend(loc='lower left', frameon=False, fontsize=7)
tag(ax, '(a)')

ax = axes[1]
groups = [
    [s for s in net['scan'] if s['label'] in ('mu0.5', 'mu1', 'mu2')],
    [s for s in net['scan'] if s['label'] in ('a0.5', 'a3')],
    [s for s in net['scan'] if s['label'] == 'N5e3'],
]
ticklabels, vals, errs, cols = [], [], [], []
for grp, col in zip(groups, [OI['blue'], OI['verm'], OI['green']]):
    for s in grp:
        vals.append(s['dev_mean']); errs.append(1.96 * s['dev_sd']); cols.append(col)
        lab = {'mu0.5': r'$\mu=0.5$', 'mu1': r'$\mu=1$', 'mu2': r'$\mu=2$',
               'a0.5': r'$\alpha=0.5$', 'a3': r'$\alpha=3$',
               'N5e3': r'$N=5{\times}10^3$'}[s['label']]
        ticklabels.append(lab)
ax.bar(range(len(vals)), vals, yerr=errs, capsize=2.5, color=cols, width=0.6)
for x in (2.5, 5.5):
    ax.axvline(x, color='0.8', lw=0.7)
ax.set_xticks(range(len(vals))); ax.set_xticklabels(ticklabels, fontsize=7.5)
ax.set_ylim(0, 19)
ax.set_ylabel(r'ledger deviation at $\langle k\rangle\simeq 8$ (%)')
tag(ax, '(b)')
fig.tight_layout()
fig.savefig(f'{FIGD}/fig7-networkci.pdf', bbox_inches='tight')
print('fig7-networkci saved')

# ================= 图级数据渲染（v0.1.2 新增, v4 稿图号） =================
def fig3_law1():
    d = _load('w3_law1_data.json')
    if d is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    cols = {0.5: OI['blue'], 1.0: OI['black'], 2.0: OI['verm']}
    for c in d['panels']:
        mu = c['mu']
        lab = fr'$\mu={mu:g}$, plateau ${c["plateau_mean"]:.3f}$'
        ax.plot(c['t'], c['scaled'], lw=0.9, color=cols[mu], label=lab)
    ax.axhline(1.0, ls='--', lw=0.7, color='0.6')
    ax.set_xscale('log'); ax.set_ylim(0.8, 1.2)
    ax.set_xlabel('$t$ (predicted tail factors applied)')
    ax.set_ylabel(r'scaled $i(t)$')
    ax.legend(loc='lower right', frameon=False, fontsize=6.5)
    tag(ax, '(a)')
    ax = axes[1]
    for c in d['panels']:
        mu = c['mu']
        ax.plot(c['i'], c['S_over_i'], lw=0.9, color=cols[mu])
    ii = np.geomspace(1e-8, 1e-1, 40)
    ax.plot(ii, 1.9 * ii ** -0.5, ':', lw=0.8, color='0.5')
    ax.plot(ii, np.full_like(ii, 1.0), ':', lw=0.8, color='0.5')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_ylim(0.5, 1e3)
    ax.set_xlabel('$i$'); ax.set_ylabel('$S/i$')
    tag(ax, '(b)')
    fig.tight_layout()
    fig.savefig(f'{FIGD}/fig3-law1.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig3 saved')


def fig4_law2():
    d = _load('w3_law2_data.json')
    if d is None:
        return
    pts = d['points']
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    eps_mk = {60: 's', 600: 'o', 6000: '^'}
    seen = set()
    for c in pts:
        g = c['alpha'] / (c['alpha'] + c['mu'])
        m = eps_mk[c['eps']]
        lab = None
        if m not in seen:
            lab = fr'$\epsilon={c["eps"]}$'
            seen.add(m)
        ax.plot(g, c['x1_mean'], m, ms=4, color=OI['blue'], label=lab)
    ax.plot([0, 1], [1, 0], '--', lw=0.8, color='0.6')
    ax.set_xlabel(r'$\alpha/(\alpha+\mu)$')
    ax.set_ylabel(r'$x_1/N$ (SSA)')
    ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    tag(ax, '(a)')
    ax = axes[1]
    seen = set()
    for c in pts:
        g = c['share2_law']
        m = eps_mk[c['eps']]
        lab = None
        if m not in seen:
            lab = fr'$\epsilon={c["eps"]}$'
            seen.add(m)
        ax.plot(g, c['share2'], m, ms=4, color=OI['verm'], label=lab)
    ax.plot([0, 1], [0, 1], '--', lw=0.8, color='0.6')
    ax.set_xlabel(r'law $\alpha/(\alpha+\mu)$')
    ax.set_ylabel(r'$x_2/(x_1+x_2)$ (SSA)')
    ax.legend(frameon=False, fontsize=6.5, loc='upper left')
    tag(ax, '(b)')
    fig.tight_layout()
    fig.savefig(f'{FIGD}/fig4-law2.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig4 saved')


def fig5_networkledger():
    d = _load('w6_network_ledger_data.json')
    if d is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    t = np.array(d['traj']['a0.5'][0]['t'])
    for key, col, lab in [('a0.5', OI['blue'], r'$\alpha=0.5$'),
                          ('a3', OI['verm'], r'$\alpha=3$')]:
        arr = np.array([r['S'] for r in d['traj'][key]])
        ax.plot(t, arr.mean(0), lw=1, color=col, label=lab)
        ax.fill_between(t, np.percentile(arr, 2.5, axis=0),
                        np.percentile(arr, 97.5, axis=0), color=col, alpha=0.2)
    ax.set_xscale('log')
    ax.set_xlabel('$t$'); ax.set_ylabel('$S(t)/N$')
    ov = d['summary']['ci_overlap_frac']
    ax.text(0.97, 0.06, f'95% CI overlap: {ov*100:.1f}%',
            transform=ax.transAxes, ha='right', fontsize=7)
    ax.legend(loc='upper right', frameon=False, fontsize=7)
    tag(ax, '(a)')
    ax = axes[1]
    summ = d['summary']
    for key, col in [('a0.5', OI['blue']), ('a3', OI['verm'])]:
        s = summ[key]
        law = s['law_share']
        obs = s['x2_share_mean']
        ax.plot([law], [obs], 'o', ms=5, color=col)
        ax.annotate(fr"$\alpha={key[1:]}$: {obs:.3f}", (law, obs),
                    textcoords='offset points', xytext=(6, -2), fontsize=7)
    ax.plot([0, 1], [0, 1], '--', lw=0.8, color='0.6')
    ax.set_xlabel(r'law $\alpha/(\alpha+\mu)$')
    ax.set_ylabel(r'network $x_2/(x_1+x_2)$')
    tag(ax, '(b)')
    fig.tight_layout()
    fig.savefig(f'{FIGD}/fig5-networkledger.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig5 saved')


def fig6_boundary():
    d = _load('w6_failure_boundary_data.json')
    if d is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    ks = [r['k'] for r in d['reach']]
    rf = [r['rf_N_mean'] for r in d['reach']]
    er = [r['rf_N_sd'] for r in d['reach']]
    hmf = [r['hmf_reach'] for r in d['reach']]
    ax.errorbar(ks, rf, yerr=er, fmt='o-', ms=4, lw=1, color=OI['blue'],
                capsize=2, label='quenched (8 seeds)')
    ax.plot(ks, hmf, '^-', ms=4, lw=1, color=OI['orange'], label='HMF')
    ax.axhline(0.999, ls='--', lw=0.7, color='0.6')
    ax.text(ks[-1], 0.969, 'MF: $R_f=1-\\epsilon$',
            ha='right', fontsize=7, color='0.4')
    ax.set_xscale('log')
    ax.set_xlabel(r'mean degree $\langle k\rangle$')
    ax.set_ylabel('$R_f/N$')
    ax.legend(loc='lower right', frameon=False, fontsize=7)
    tag(ax, '(a)')
    ax = axes[1]
    ax.errorbar(ks, [r['dev_mean'] for r in d['ledger_dev_poisson']],
                yerr=[r['dev_sd'] for r in d['ledger_dev_poisson']],
                fmt='o-', ms=4, lw=1, color=OI['blue'], capsize=2,
                label='Poisson')
    for r in d['ledger_dev_sf']:
        ax.errorbar([r['k_matched']], [r['dev_mean']], yerr=[r['dev_sd']],
                    fmt='D', ms=5, color=OI['verm'], capsize=2,
                    label=fr'scale-free $\tau={r["tau"]:g}$')
    xx = np.array(ks, float)
    sl = d['meta'].get('fit_exponent', -0.71)
    ax.plot(xx, np.array([r['dev_mean'] for r in d['ledger_dev_poisson']])[0] *
            (xx / xx[0]) ** sl, '--', lw=0.8, color='0.5')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'mean degree $\langle k\rangle$')
    ax.set_ylabel('ledger deviation (%)')
    ax.legend(loc='upper right', frameon=False, fontsize=6.5)
    tag(ax, '(b)')
    fig.tight_layout()
    fig.savefig(f'{FIGD}/fig6-boundary.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig6 saved')


def fig8_lonespreader():
    d = _load('w4_lone_spreader_data.json')
    if d is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    Ns = sorted({c['N'] for c in d['cells']})
    mus = sorted({c['mu'] for c in d['cells']})
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(mus)))
    for mu, col in zip(mus, cols):
        ys = [next(c['upper_bound_95'] for c in d['cells']
                   if c['N'] == N and c['mu'] == mu) for N in Ns]
        ax.plot(Ns, ys, 'o-', ms=3.5, lw=1, color=col, label=fr'$\mu={mu:g}$')
    ax.plot(Ns, [3 / n for n in Ns], '--', lw=0.8, color='0.6')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('$N$'); ax.set_ylabel(r'95% upper bound on $P(i_f>0)$')
    ax.legend(frameon=False, fontsize=6, ncol=2, loc='lower left')
    tag(ax, '(a)')
    ax = axes[1]
    for mu, col in zip(mus[:1] + mus[-1:], [cols[0], cols[-1]]):
        xs = [1 / N for N in Ns]
        ys = [next(c['one_minus_Rf'] for c in d['cells']
                   if c['N'] == N and c['mu'] == mu) for N in Ns]
        ax.plot(xs, ys, 'o-', ms=3.5, lw=1, color=col, label=fr'$\mu={mu:g}$')
    ax.plot([1 / max(Ns), 1 / min(Ns)], [1 / max(Ns), 1 / min(Ns)],
            '--', lw=0.8, color='0.6')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('$1/N$'); ax.set_ylabel(r'ledger correction $1-R_f/N$')
    ax.legend(frameon=False, fontsize=7, loc='upper left')
    tag(ax, '(b)')
    fig.tight_layout()
    fig.savefig(f'{FIGD}/fig8-lonespreader.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig8 saved')


def fig9_realcascades():
    m = _load('w5_main_results.json')
    if m is None:
        return
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.0))
    ax = axes[0]
    for blk, col in [(b, OI['blue']) for b in m['theta']]:
        th = blk['theta_h']
        ax.errorbar([th], [blk['ic_med']],
                    yerr=[[blk['ic_med'] - blk['ic_ci'][0]],
                          [blk['ic_ci'][1] - blk['ic_med']]],
                    fmt='o', ms=4, color=col, capsize=2)
    ax.axhline(1.0, ls='--', lw=0.8, color='0.6')
    ax.set_xlabel(r'session threshold $\theta$ (h)')
    ax.set_ylabel('regression intercept $a$')
    ax.set_ylim(0, 1.6)
    tag(ax, '(a)')
    ax = axes[1]
    rep = None
    for r in m['per_cascade_RT'].get('theta_6h', []):
        if abs(r.get('n_users', 0) - 661) <= 30 and r.get('r2', 0) > 0.5:
            rep = r
            break
    if rep is not None:
        ax.plot(rep['px'], rep['py'], 'o', ms=3, color=OI['blue'])
        xx = np.array(rep['px'], float)
        ax.plot(xx, rep['intercept'] + rep['mu'] * xx, '-', lw=1,
                color=OI['blue'],
                label=fr'fit: $a={rep["intercept"]:.2f}$, $\mu={rep["mu"]:.2f}$')
        ax.plot(xx, 1 - rep['mu'] * xx, '--', lw=0.8, color='0.6',
                label='mean field')
        ax.legend(frameon=False, fontsize=7)
        ax.set_xlabel('$S/i$'); ax.set_ylabel(r'$\Delta S/\Delta i$')
        tag(ax, f'(b) $P={rep["n_users"]}$, $R^2={rep["r2"]:.2f}$')
    ax = axes[2]
    for blk, col in [(b, OI['blue']) for b in m['theta']]:
        th = blk['theta_h']
        lo = blk['alpha_pooled_ci'][0]
        hi = blk['alpha_pooled_ci'][1]
        ax.plot([th, th], [lo, hi], '-', lw=2.4, color=col)
    re = m['re_view']
    ax.plot([8, 8], re['alpha_pooled_ci'], '-', lw=2.4, color=OI['verm'],
            label='RE view')
    ax.set_yscale('log')
    ax.set_xlabel(r'session threshold $\theta$ (h; RE at 8)')
    ax.set_ylabel(r'$\hat\alpha$ pooled, cluster bootstrap 95% CI')
    tag(ax, '(c)')
    fig.tight_layout()
    fig.savefig(f'{FIGD}/fig9-realcascades.pdf', bbox_inches='tight')
    plt.close(fig)
    print('fig9 saved')


for _f in (fig3_law1, fig4_law2, fig5_networkledger, fig6_boundary,
           fig8_lonespreader, fig9_realcascades):
    _f()
print('all extra figures done')
