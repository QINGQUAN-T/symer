# -*- coding: utf-8 -*-
"""W7-R5 测量模型数据包（审稿意见 #2；v0.1.2 重写为自含管线）。

对 Higgs-Twitter RT/RE 档案、theta∈{2,6,24}h 输出:
  - 守卫漏斗: >=30事件源数 -> span>=2theta -> 休眠删失守卫 -> 回归过箱
  - 原始计数: 单爆发/多爆发/账本池; 每级联末位传播者删失数(=级联数)
  - x2_pooled 及其簇 bootstrap CI
  - mu_hat / identity截距 / R2 中位数, mu_ci / R2 IQR / mu IQR（每级联分布）
  - 每级联 alpha_hat 分布（中位数、四分位、bootstrap CI）
  - alpha_pooled = mu_med * x2/x1 及其联合簇 bootstrap CI（审稿意见 M1:
    池化估计量与每级联估计量分开报告，不确定性由级联级联合重采样给出）
  - RE(回复)视角: theta=6h 同一管线（论文 §7.2(iii)）

输出: results/w7_measure.json（汇总）
      results/w5_main_results.json（逐级联记录 + RT/RE 汇总, 供 w7_figs 图 9）
随机种子: bootstrap 固定 seed=20261018, B=2000。
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w5_invert import load_rt, analyze_cascade

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/results/w7_measure.json'
OUT_MAIN = f'{ROOT}/results/w5_main_results.json'


def find_data(rtype='RT'):
    """按 data_README.md 约定寻找数据文件。"""
    for name in ('higgs-activity.txt.gz', 'higgs-activity_time.txt.gz'):
        p = f'{ROOT}/data/{name}'
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        'Higgs activity_time 数据未找到。请按 data_README.md 下载到 data/ 目录。')


def guarded_cascades(casc, th, t_end):
    """与 W5/W7 完全相同的三级守卫; 返回 (rows, funnel)。"""
    n0 = sum(1 for evs in casc.values() if len(evs) >= 30)
    funnel = dict(n_sources_ge30=n0, dropped_span=0, dropped_dormant=0,
                  passed_guards=0, passed_regression=0)
    rows = []
    for v, evs in casc.items():
        if len(evs) < 30:
            continue
        ts = [t for _, t in evs]
        if max(ts) - min(ts) < 2 * th:
            funnel['dropped_span'] += 1
            continue
        if max(ts) > t_end - 24 * 3600.0:
            funnel['dropped_dormant'] += 1
            continue
        funnel['passed_guards'] += 1
        c = analyze_cascade(evs, th, 16)
        if c is not None:
            funnel['passed_regression'] += 1
            c['n_events'] = len(evs)
            rows.append(c)
    return rows, funnel


def cluster_bootstrap(rows, seed=20261018, boot=2000):
    """级联级联合簇 bootstrap（审稿意见 M1）。

    每次重采样: 有放回抽级联 -> mu_med = median(mu_c);
    x1/x2 池化计数求和 -> alpha_pool = mu_med * (sum x2)/(sum x1)。
    对每级联中位数与池化估计量给出联合置信区间。"""
    rng = np.random.default_rng(seed)
    n = len(rows)
    mus = np.array([r['mu'] for r in rows])
    # 池化计数（注意: analyze_cascade 的 x1/x2 字段是每级联"分数"，
    # 池化估计必须用 n_bursts 计数字段, 与 w7_negctl/w5 演示的口径一致）
    x1 = np.array([r['n_bursts'][1] for r in rows], float)
    x2 = np.array([r['n_bursts'][2] for r in rows], float)
    als = np.array([r['alpha'] for r in rows])
    x2s = np.array([r['x2'] for r in rows], float)
    ics = np.array([r['intercept'] for r in rows])
    mu_bs, pool_bs, x2p_bs, med_al, med_x2, med_ic = [], [], [], [], [], []
    for _ in range(boot):
        idx = rng.integers(0, n, n)
        m_b = np.median(mus[idx])
        x1_b, x2_b = x1[idx].sum(), x2[idx].sum()
        mu_bs.append(m_b)
        pool_bs.append(m_b * x2_b / x1_b)
        x2p_bs.append(x2_b / (x1_b + x2_b))   # 池化份额 = multi/(single+multi)
        med_al.append(np.median(als[idx]))        # 每级联 alpha 中位数
        med_x2.append(np.median(x2s[idx]))
        med_ic.append(np.median(ics[idx]))
    q = lambda a, p: float(np.percentile(a, p))
    return dict(
        mu_med=float(np.median(mus)),
        mu_ci=[q(mu_bs, 2.5), q(mu_bs, 97.5)],
        mu_q1=float(np.percentile(mus, 25)), mu_q3=float(np.percentile(mus, 75)),
        ic_med=float(np.median(ics)),
        ic_ci=[q(med_ic, 2.5), q(med_ic, 97.5)],
        r2_med=float(np.median([r['r2'] for r in rows])),
        r2_q1=float(np.percentile([r['r2'] for r in rows], 25)),
        r2_q3=float(np.percentile([r['r2'] for r in rows], 75)),
        x2_pooled=float(x2.sum() / (x1.sum() + x2.sum())),
        x2_pooled_ci=[q(x2p_bs, 2.5), q(x2p_bs, 97.5)],
        x2_med=float(np.median(x2s)),
        x2_med_ci=[q(med_x2, 2.5), q(med_x2, 97.5)],
        alpha_percasc=dict(
            median=float(np.median(als)),
            q1=float(np.percentile(als, 25)), q3=float(np.percentile(als, 75)),
            ci=[q(med_al, 2.5), q(med_al, 97.5)]),
        alpha_pooled=float(np.median(mus) * x2.sum() / x1.sum()),
        alpha_pooled_ci=[q(pool_bs, 2.5), q(pool_bs, 97.5)],
        boot=boot, seed=seed)


def summarize_rt(casc, t_end, thetas=(2.0, 6.0, 24.0)):
    res, per = {}, {}
    for th_h in thetas:
        th = th_h * 3600.0
        rows, funnel = guarded_cascades(casc, th, t_end)
        boot = cluster_bootstrap(rows)
        single = sum(r['n_bursts'][1] for r in rows)
        multi = sum(r['n_bursts'][2] for r in rows)
        block = {
            'dropped_span': funnel['dropped_span'],
            'dropped_dormant': funnel['dropped_dormant'],
            'passed_guards': funnel['passed_guards'],
            'passed_regression': funnel['passed_regression'],
            'gate_pass_rate': round(funnel['passed_regression'] /
                                    max(1, funnel['passed_guards']), 3),
            'raw_single_burst': single,
            'raw_multi_burst': multi,
            'ledger_pool': single + multi,
            'lone_spreader_censored': funnel['passed_regression'],
            'x2_pooled': boot['x2_pooled'],
            'x2_pooled_ci': boot['x2_pooled_ci'],
            'mu_med': boot['mu_med'],
            'mu_ci': boot['mu_ci'],
            'mu_q1': boot['mu_q1'], 'mu_q3': boot['mu_q3'],
            'ic_med': boot['ic_med'],
            'ic_ci': boot['ic_ci'],
            'r2_med': boot['r2_med'],
            'r2_q1': boot['r2_q1'], 'r2_q3': boot['r2_q3'],
            'x2_med': boot['x2_med'], 'x2_med_ci': boot['x2_med_ci'],
            'alpha_percasc': boot['alpha_percasc'],
            'alpha_pooled': boot['alpha_pooled'],
            'alpha_pooled_ci': boot['alpha_pooled_ci'],
            'n_users_med': int(np.median([r['n_users'] for r in rows])),
            'events_med': int(np.median([r['n_events'] for r in rows])),
            'boot': boot['boot'], 'seed': boot['seed'],
        }
        key = f'theta_{th_h:g}h'
        res[key] = block
        per[key] = rows
        print(f"  theta={th_h}h: n={block['passed_regression']} "
              f"x2={boot['x2_pooled']:.4f} mu={boot['mu_med']:.3f} "
              f"ic={boot['ic_med']:.2f} a_pool={boot['alpha_pooled']:.4f} "
              f"CI=[{boot['alpha_pooled_ci'][0]:.4f},{boot['alpha_pooled_ci'][1]:.4f}]")
    return res, per


def main():
    data_path = find_data('RT')
    print('data:', os.path.basename(data_path))
    casc = load_rt(data_path, 'RT')
    t_end = max(t for evs in casc.values() for _, t in evs)
    print('RT events:', sum(len(v) for v in casc.values()),
          '| sources:', len(casc), '| t_end(h):', round(t_end / 3600, 2))

    res, per = summarize_rt(casc, t_end)
    json.dump(res, open(OUT, 'w'), indent=1, ensure_ascii=False)
    print('->', os.path.relpath(OUT, ROOT))

    # RE(回复)视角: theta=6h（论文 §7.2(iii)，n=8 探索性）
    casc_re = load_rt(data_path, 'RE')
    t_end_re = max(t for evs in casc_re.values() for _, t in evs)
    print('RE events:', sum(len(v) for v in casc_re.values()),
          '| sources:', len(casc_re))
    res_re, per_re = summarize_rt(casc_re, t_end_re, thetas=(6.0,))

    main_out = {
        'meta': {
            'dataset': 'SNAP Higgs-Twitter activity_time',
            'pipeline': 'w5_invert.analyze_cascade (bins=16), guards: >=30 events, '
                        'span>=2theta, dormancy censor 24h; ledger pool = P-1 (Theorem 8)',
            'bootstrap': 'cascade-level cluster bootstrap, B=2000, seed=20261018',
            'version': 'symer v0.1.2',
        },
        'theta': [{'theta_h': th_h, **res[f'theta_{th_h:g}h']} for th_h in (2.0, 6.0, 24.0)],
        're_view': {'theta_h': 6.0, **res_re['theta_6h']},
        'per_cascade_RT': {k: [{kk: vv for kk, vv in r.items()} for r in v]
                           for k, v in per.items()},
        'per_cascade_RE_6h': per_re['theta_6h'],
    }
    json.dump(main_out, open(OUT_MAIN, 'w'), indent=1, ensure_ascii=False)
    print('->', os.path.relpath(OUT_MAIN, ROOT))


if __name__ == '__main__':
    main()
