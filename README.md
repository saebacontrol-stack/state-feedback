# State-Feedback for 2 mass system Demo

Pythonのcontrolライブラリおよびnumpy、matplotlib を用いて、2慣性系モータモデルに対する連続時間設計とディジタル（離散時間）設計の状態フィードバック制御を比較するためのシミュレーションスクリプトです。

## 概要 (Overview)
モータと負荷がバネ（減速機剛性）で結合された「2慣性系モデル」を対象に、状態フィードバック制御をかけたときのモータ速度の応答を比較します。  
ディジタル制御のサンプリング周期（ $T_s$）と、微小時間ステップ（ $dt$、4次ルンゲ・クッタ法による数値積分）を分けることで、実際のマイコン制御に近い挙動を再現しています。

---

## 実行方法 (Usage)

### 必要なライブラリ
```bash
pip install numpy matplotlib control
```
### 実行手順
```bash
python sample_state-feedback.py
```
