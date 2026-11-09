# electron・muonの温度別run統合

2023年PS/SPSの保存済みADC分布を、同じbeam・粒子・channel・温度区間でrunを跨いで加算するC++/ROOTコードです。PSとSPS、electronとmuonは必ず別にします。ユーザー確認に基づき、bias電圧とASIC gain設定はPS内・SPS内でそれぞれ共通と宣言して使用します。

元のrawを読み直さず、`analysis/result/sipm_gain/<from-label>/data/{electron,muon}/<sample>/collection/spectra.root` の `Spectra` とADC histogramを読みます。元の0.5℃区間を維持し、memory cell 0・HitTag 1・GainTag 1・温度情報ありの分布を使用します。2 run以上が寄与する分布だけをfit用へ渡します。温度幅が異なる入力は拒否します。温度の平均・RMSは各入力の有効温度entry数で統合します。

同じ粒子内では異なるbeam energyのrunも含めます。元run・energy・元sid・寄与entry数は `Contributions` treeに保存します。元のADC軸で加算し、run間のpedestal位置合わせは行いません。全体分布の形やピークの位置がrun間で異なる場合、統合でピークが鈍ることもあるため、統合後のfit品質・4手法の一致を確認してください。統合だけでgainが必ず決まるわけではありません。

## 保存先

- `analysis/result/sipm_gain_pool/<label>/inputs/`：入力ROOTのchecksum・fingerprint、設定、ソース控え。
- `analysis/result/sipm_gain_pool/<label>/data/<ps|sps>/<electron|muon>/layer_XX/`：統合ADC、元runの寄与、温度を℃で記した一覧、加算数の検証記録。
- `analysis/result/sipm_gain/<label>/data/<particle>/<beam>_<particle>_temperature_pool/collection/`：既存fitコードが読める最終ADC ROOT。
- 以後の手法別gain・比較ROOT・図は既存の `analysis/result/sipm_gain/<label>/` 配下。元labelの結果は変更しません。
- pooling・compare・plotsのSlurmログは `analysis/log/sipm_gain_pool/<label>/<stage>/`。既存fitのログは `analysis/log/sipm_gain/<label>/fit/`。

各分布に寄与したrun数を `Spectra.pooled_runs` に保存し、最終ROOTにも `Contributions` を保持します。異なる温度区間やPS/SPSを加算しないこと、ADC countの保存、入力変更を検査します。prepareで入力ROOTをchecksum検証し、collectではそのfingerprintと完了印・収集設定のchecksumを再確認します。

## 実行方法

ROOT環境でcompileします。既存 `sipm_gain/Analyze` も必要です。

```bash
cd /home/murata_t/ScECAL_BeamTest
make -C analysis/sipm_gain_pool -j2
```

以下は**1段階ずつ、前のジョブの正常終了を確認してから**実行します。1コマンドが次の段階まで自動投入することはありません。既定の新しいlabelは `full_grouped_temperature_pool`、入力は `full_grouped` です。

```bash
# 1. 92 runの入力検証・設定の固定（1 job）
bash analysis/sipm_gain_pool/submit_stage.sh prepare

# 2. PS/SPS × electron/muon × 30 layerでADCを統合（120 task、最大32並列）
bash analysis/sipm_gain_pool/submit_stage.sh collect

# 3. 30 layer分を各beam・粒子のfit用ROOTへまとめる（4 task）
bash analysis/sipm_gain_pool/submit_stage.sh merge

# 4. 統合ADCへ既存の4手法を適用（480 task、最大32並列）
bash analysis/sipm_gain_pool/submit_stage.sh fit

# 5. 4手法の数値比較（1 job）
bash analysis/sipm_gain_pool/submit_stage.sh compare

# 6. 比較plot（30 task、最大32並列）
bash analysis/sipm_gain_pool/submit_stage.sh plots
```

labelを変える場合は全段階で同じ第2引数を指定します。第3引数で入力labelを指定できます。`submit_stage.sh` は全PS/SPS・全粒子・全30 layer用です。部分検証では下の直接コマンドを使ってください。条件やコードを変える場合は、新しいlabelにします。

## 小規模な直接実行

```bash
./analysis/sipm_gain_pool/Pool --from-label full_grouped --label temperature_pool_validation --stage prepare --settings-policy common-within-beam --beam sps --kind muon --layers 0 --channels 20,21,20010,20011
./analysis/sipm_gain_pool/Pool --from-label full_grouped --label temperature_pool_validation --stage collect --beam sps --kind muon --layer 0
./analysis/sipm_gain_pool/Pool --from-label full_grouped --label temperature_pool_validation --stage merge --beam sps --kind muon
./analysis/sipm_gain/Analyze run --manifest analysis/result/sipm_gain/temperature_pool_validation/inputs/manifest.tsv --label temperature_pool_validation --stage fit --channel-groups 1 --plot-limit 4
```

この小規模labelは全channelのcalibrationではありません。pooling後の結果は、まず手法別gainの候補と診断です。既存LEDの採用結果へ自動的に混ぜたり、calibrationの未決定channelへ無条件に代入したりしません。温度統合後の手法一致と元run間の安定性を確認してから、採用条件を決めるための入力になります。LEDの2DAC条件による採用ROOTは、独立した `sipm_gain_selection` の段階で作成します。
