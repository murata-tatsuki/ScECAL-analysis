# hl_intercept: HG–LG intercept修正の全channel比較

`channel_response/adc_energy` が保存したROOTを読み、data・修正前MC・修正後MCを同じHG範囲で比較する。元のdecode/calib ROOT、本番の校正定数、既存の集計ROOTは変更しない。既定は100 GeV、接続済み6300 physical channel（layer 0–29、chip 5はchannel 0–29、他は0–35）。bad-channelも診断対象に含める。

## 実行

```bash
cd ~/ScECAL_BeamTest/analysis/channel_response/hl_intercept
make
make test
bash run.sh
```

入力は `../../result/channel_response/hl_intercept/` の次の4ファイル。

- `after_correction/data/100GeV.root`
- `before_correction/data/100GeV.root`（dataが同じであることの検証用）
- `before_correction/simulation/threshold/100GeV.root`
- `after_correction/simulation/threshold/100GeV.root`

`--base DIRECTORY`、`--reference ROOT`、`--energy 100`で変更できる。新しいenergyを指定した場合、そのenergyの4入力が必要。energy間の合算はしない。

```bash
bash run.sh --help
./CompareHLIntercept --help
```

既存出力を上書きしない。再計算や異なるfit範囲の比較は、別のbaseに入力2ディレクトリへのsymlinkを用意して実行する。`--channels 240008,420008`は少数channelの開発確認用で、既定は全channel。失敗時は隠しstagingディレクトリを残し、その場所を表示する。正常終了時のみ完成した4出力を保存先へ移動する。

## 出力

`analysis/result/channel_response/hl_intercept/` 直下:

```text
InterCalib_data_100GeV.root
InterCalib_before_100GeV.root
InterCalib_after_100GeV.root
comparison/
  100GeV/
    channels.tsv           # 全channel × 3sample、fit・残差・hit・energy診断
    comparison.tsv         # 全channel、dataとの差と修正前後差
    summary.txt            # fit可否、改善/悪化channel数、集計
    summary.root           # map、係数差のTGraphErrors、分布、overview canvas
    overview.png
    overview.pdf
    channels/cell_*.png     # 不一致が残る/悪化が大きいchannelの個別図
    provenance.txt         # 入力、event数、fit条件
    command.txt
    checksums.txt
    input_file_stats.txt
    run.log
```

個別図は既定で「afterの残存不一致が大きい20 channel」と「悪化量が大きい20 channel」の和集合。`--top`で変更できる。数値とmapは個別図の対象にかかわらず全channelを含む。mapの未評価セルはNaNで、`valid_common_fit`も別途保存する。改善/悪化の区分は差の符号による記述で、有意差判定ではない。

## Channel別の図

各channelの図は左上がHG–LG profileとfit、右上がprofile残差、左下がhit残差分布。右下はinterceptとslopeの2パネルを並べる。両パネルの横軸は `data` / `before` / `after`、色は黒 / 橙 / 青で、fitの形式誤差も表示する。破線は参照ROOTの入力係数。各パネル上部に `after - before` の係数変化を数値で表示する。

## InterCalib treeの互換性

参照ファイルの既定値:

`/home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/share/all_hl_electron2023_hlratio_v3.root`

3 ROOTに以下の同名・同型の`InterCalib` TTreeを保存する。実行時に参照treeのschemaを検証する。

| branch | ROOT leaf型 | 今回の内容 |
|---|---|---|
| CellID | `/I` | physical channel ID |
| Slope | `/D` | pedestal差引き後 `LG = Slope * HG + Intercept` |
| Intercept | `/D` | LG ADC単位の切片 |
| SlopeError | `/D` | profileの標準誤差を用いたWLSの形式誤差 |
| InterceptError | `/D` | 同上 |
| ChiSquare | `/D` | 使用binに対するWLSのchi-square |
| NDF | `/D` | 使用bin数 − 2（参照と同じdouble） |
| XMax | `/D` | **参照ROOTの値を継承** |
| Statistics | `/I` | 共通fit bin内のhit数 |

`XMax`は元のMC digitizationでHG切替に使われる値であり、今回測定した非線形開始位置でもfit上限でもない。参照コードのmultigaus fitを再現したものではなく、**tree schemaの互換性**を維持した新しいprofile fitである。3 ROOTは診断出力で、本番校正への自動適用はしない。

`InterCalib`には数値的に有効なfitのみ保存する。統計不足・欠落channelを偽の係数や入力係数で埋めない。全channelの可否は別tree `FitDiagnostics` とTSVに保存する。chi2/NDFが大きいfitは削除せず、値を保存しsummaryでも数える。校正値として採用する際はfit範囲依存性・非線形性・裾の影響を別に評価する。

`FitDiagnostics`にはstatus、bin数、hit数、実際のHG bin中心の最小/最大、slope–intercept共分散、基準HGでの予測LGと誤差を保存する。基準HGが使用範囲外なら外挿なので、範囲情報と合わせて読む。

Status: `0=valid`, `1=missing_channel`（3sampleに共通のchannelがない場合も含む）, `2=insufficient_common_bins`, `3=insufficient_HG_span`, `4=insufficient_hits`, `5=invalid_fit`, `6=missing_reference`。

## Fitの方法と既定条件

- `hg_lg_zoom_profile`を用い、bin中心HGに対する平均LGを誤差の逆二乗で重み付けした直線fit。平均LGにはTH2のY表示範囲外のhitも含まれる。
- data/before/afterの**3sample共通bin**を使う。binの全幅がfit範囲内にあり、各sampleで20 hit以上、有限で正の標準誤差があるbinのみ。
- 要求範囲はHG pedestal差引き後800–2200 ADC。MCのraw HG切替が `reference.XMax − 600` である実装に合わせ、pedestalと1 HG bin幅を引いた位置より下へ上限を制限する。`--switch-margin`はこの600を変更する。生成時の設定が異なる場合は明示的に変更する。
- 8 bin以上、bin中心の幅400 ADC以上、各sampleの採用hit合計200以上を要求。fit値が非有限、slopeが非正、または参照XMaxがない場合は校正treeに出さない。
- slopeとinterceptの誤差・共分散はWLSの形式誤差。chi2で水増ししない。binning、fit範囲、HGの測定誤差、外れ値由来の系統誤差は含まれない。
- 内側の1000–2000 ADCでも同じ手順でfitし、係数と基準HG=1500 ADCでの予測値の変化をTSVへ保存。内側fitの統計不足も別statusで示す。
- TProfileにはHG binごとの厳密な平均HGがないので、HGはbin中心近似。精密なintercalibrationには元rawを再読込し、HGの分布・誤差と裾を扱うfitが必要。

## 比較指標

1. **係数の変化と入力への閉じ**: `after.slope − before.slope`、`after.intercept − before.intercept − input.intercept`、`before.intercept`、`after.intercept − input.intercept`、`after.slope − input.slope`。期待は直線領域で概ね0。dataとの差と基準HGでの差も保存。
2. **fit非依存のdata一致度**: 共通HG binで `MC平均LG − data平均LG` のbiasとRMSを計算。bin重みは `min(Ndata,Nbefore,Nafter)` でbefore/afterに同じものを使う。RMS改善量は `RMS_before − RMS_after`（正が改善）。共通binがあればfit失敗channelにも指標を残す。summaryの集約は3fit有効channelのみ。
3. **残差の幅と裾**: TH2のbin中心から `LG − fitted LG` を作り、平均、RMS、中央値、中央68%の半幅を計算。各sampleの中央値からの距離がdataの68%半幅の5倍を超える割合を同じ閾値で比較する。これは有限bin幅による近似で、Y under/overflowは除外し、その数を別列で示す。profile fitはY flowを含むため、両者の母集団は厳密には同じでない。
4. **energyへの影響**: `hg_selected_energy_profile`で同じHG領域の最終採用energyをdataと比較。TH2内の平均energy、energy和/event、表示範囲外率も出力。範囲外hitのenergy値は復元できないため、後者は必ず`inrange`と表記する。
5. **HG/LG接続の診断**: raw HG=2600（`--gain-switch`で変更）から各sampleのpedestalを引き、その左右250 ADCで局所直線fit。境界での右−左予測energyを`switch_jump_proxy`として保存する。左右各3 bin・60 ADC幅・60 hit以上が必要。MCのHG置換で右側に点がない場合は未評価。これは選択branchや同一hitのHG/LG差を直接測ったものではない。
6. **hitとADC分布**: hit数/event、範囲内HG/LG平均、flow割合、before/afterのHGおよびLGの正規化射影のtotal variationを出力。TVは0が同一、1が重なりなしで、軸のflowを含む。独立に生成したMCではHG分布も統計的に揺らぐため、TV=0を必須にしない。

単一channelのTH2/Profileにはevent相関が保存されていない。event全energyのresolution、正確なLG採用率、同一hitの `E_LG−E_HG` はここからは算出しない。

dataの2コピーは、入力manifest、選別metadata、event数、全6種類のTH2/Profileの全bin・誤差・統計量・pedestalが一致することを検証する。不一致なら停止する。before/after MC同士を同一eventとして扱わず、係数差の図の誤差は独立sample近似である。

## 検証

`make test`は既知のHG–LG直線を持つ合成ROOTを作り、傾きと切片の復元、入力切片との差分、参照と同じ9 branch、XMax継承、統計不足・欠落channel、dataの不一致検出、出力上書き防止を検証する。

## コードの配置

解析コードは `analysis/channel_response/hl_intercept/` に配置する。共通ROOTビルド設定は `../root.mk` を参照する。出力は `analysis/result/channel_response/hl_intercept/`。保存済みの実行ログ・checksumsには、生成時の旧コードパスが残る。再実行にはこのREADMEの現行パスを使用する。
