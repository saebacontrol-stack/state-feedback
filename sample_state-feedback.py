import control as ct
import matplotlib.pyplot as plt
import numpy as np


dt = 10e-6                  # タイムステップ 10μs
t_end = 0.01                # シミュレーション時間 0.01s
steps = int(t_end/dt)       # ステップ数

Ts = 100e-6         # 制御周期[s]
Ns = int(Ts/dt)     # ディジタル系１ステップあたりの連続時間系ステップ数


t_history = np.linspace(0, t_end, steps)    # 時間データ
x_history1 = np.zeros((steps, 3))           # 状態変数データ１[steps x 3]
x_history2 = np.zeros((steps, 3))           # 状態変数データ２[steps x 3]
x_history3 = np.zeros((steps, 3))           # 状態変数データ３[steps x 3]

Jm = 1e-4       # モータイナーシャ[kgm^2]
Jl = 10.0       # リンクイナーシャ[kgm^2]
Ng = 100        # 減速比[-]
Kc = 3.0e4      # 減速機剛性[Nm/rad]

vm0 = 0.2       # 初期モータ速度[rad/s]
vL0 = vm0/Ng    # 初期リンク速度[rad/s]
xs0 = 0.0       # 初期ねじれ角[rad]

x0_real = np.array([vm0, vL0, xs0])     # 初期状態変数（実機）
x0_obs = np.array([vm0, 0.0, 0.0])      # 初期状態変数（オブザーバ）

Ac = np.array(
        [
            [0.0,  0.0,  -Kc/(Jm*Ng)],
            [0.0,  0.0,  Kc/Jl],
            [1/Ng, -1.0, 0.0],
        ]
    )                                   # 制御対象状態方程式A行列（連続時間系）
Bc = np.array([[1/Jm], [0.0], [0.0]])   # 制御対象状態方程式B行列（連続時間系）
Cc = np.array([1.0, 0.0, 0.0])          # 制御対象状態方程式C行列（連続時間系）
Dc = 0.0                                # 制御対象状態方程式D行列（連続時間系）

sys_Ps = ct.StateSpace(Ac, Bc, Cc, Dc)                  # 制御対象状態方程式（連続時間系）
sys_Pz = ct.sample_system(sys_Ps, Ts, method='zoh')     # 制御対象状態方程式（離散時間系）

Ad = sys_Pz.A   # 離散時間状態方程式A行列
Bd = sys_Pz.B   # 離散時間状態方程式B行列
Cd = sys_Pz.C   # 離散時間状態方程式C行列
Dd = sys_Pz.D   # 離散時間状態方程式D行列


pfb_s = np.array([-1500.0, -1500.0, -1500.0])       # 状態フィードバックの極[rad/s]（連続時間系）
K_c = ct.acker(Ac, Bc, pfb_s)                       # 極配置法による状態フィードバックゲインの計算（連続時間系）
         
pfb_z = np.exp(pfb_s*Ts)            # 連続時間系の極を離散時間系の極へ変換
K_d = ct.acker(Ad, Bd, pfb_z)       # 極配置法による状態フィードバックゲインの計算（離散時間系）


def sys_calc(x, u):
    """状態方程式 dx/dt = A@x + B@u 計算"""
    return Ac @ x + (Bc @ np.array([[u]])).flatten()


def solve_rk4(sys_calc, x, u, dt):
    """4次ルンゲ・クッタ法による1ステップ数値積分
    引数:
    -----------
    sys_calc : 状態方程式 dx/dt = A@x + B@u 計算用関数
    x : 現在の状態変数
    u : 現在の入力
    dt : 積分時間刻み幅[s]

    戻り値:
    --------
    x_next : dt 秒後の状態変数
    """
    k1 = sys_calc(x, u)
    k2 = sys_calc(x + 0.5 * dt * k1, u)
    k3 = sys_calc(x + 0.5 * dt * k2, u)
    k4 = sys_calc(x + dt * k3, u)

    x_next = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return x_next


def simulation():
  """シミュレーション実行関数"""
  global x_history1, x_history2, x_history3
  
  x1 = x0_real      # 初期状態の設定（状態変数１）
  x2 = x0_real      # 初期状態の設定（状態変数２）
  x3 = x0_real      # 初期状態の設定（状態変数３）

  cnt = Ns      # カウント変数

  for k in range(steps):
    x_history1[k] = x1
    x_history2[k] = x2
    x_history3[k] = x3

    u1 = -K_c@x1.T      # 状態フィードバックによる制御入力決定（連続時間系）

    if cnt==Ns:
       # 状態フィードバックによる制御入力決定（離散時間系、Nsステップで1度の制御入力更新により離散時間シミュレーションを実現）
       u2 = -K_d@x2.T       # ディジタル設計ゲイン使用
       u3 = -K_c@x3.T       # アナログ設計ゲイン使用
       cnt = 0              # カウント変数の初期化

    x1 = solve_rk4(sys_calc, x1, u1, dt)    # 1ステップ更新（制御対象：アナログ、制御器：アナログ設計）
    x2 = solve_rk4(sys_calc, x2, u2, dt)    # 1ステップ更新（制御対象：ディジタル、制御器：ディジタル設計）
    x3 = solve_rk4(sys_calc, x3, u3, dt)    # 1ステップ更新（制御対象：ディジタル、制御器：アナログ設計）

    cnt += 1    # カウントアップ


# グラフ描画
simulation()
plt.plot(t_history, x_history1[:, 0], label="Plant:analog, Controller:analog")
plt.plot(t_history, x_history2[:, 0], label="Plant:digital, Controller:digital")
plt.plot(t_history, x_history3[:, 0], label="Plant:digital, Controller:analog")
plt.xlabel("Time [s]")
plt.ylabel("Motor Velocity[rad/s]")
plt.grid(True)
plt.legend()
plt.show()