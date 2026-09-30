# Channel response解析

解析コード、ジョブログ、結果を同じ解析名のsubdirectoryにまとめています。

```text
analysis/channel_response/
  execute.sl              両解析共通のSlurm実行wrapper
  root.mk                 共通のC++・ROOTビルド設定
  adc_energy/             HG–LG・ADC–energy比較のC++・shell
  hg_lg_tail/             tailとevent内同時発生のC++・shell
  job/
    adc_energy/           応答比較のジョブログ・投入記録
    hg_lg_tail/           tail解析のジョブログ・投入記録

analysis/result/channel_response/
  adc_energy/             応答比較のROOT・PNG・入力設定
  hg_lg_tail/             tail解析のROOT・PNG・集計表・入力設定
```

- [adc_energyの説明・実行方法](adc_energy/README.md)
- [hg_lg_tailの説明・実行方法](hg_lg_tail/README.md)

共通のSlurm実行wrapperは直下の`execute.sl`、C++・ROOTのビルド設定は`root.mk`にまとめています。各解析の`submit.sh`・`Makefile`がこれらを参照します。解析ごとのジョブ名・ログ保存先はそれぞれの`submit.sh`で設定します。

`FastPng.hh`は`adc_energy/`だけで使用するため、その中に置いています。tailの校正係数読取は`analysis/impact_studies/StudyCommon.hh`を参照します。

## 実行例

```bash
# channel別のHG–LG・ADC–energy比較：既存ROOTから再描画
cd ~/ScECAL_BeamTest/analysis/channel_response/adc_energy
bash submit.sh --energy 100 --plots-only

# HG–LG tailとevent内同時発生
cd ~/ScECAL_BeamTest/analysis/channel_response/hg_lg_tail
bash submit.sh --energy 100
```

`adc_energy`の`--plots-only`は既存の集計ROOTを使います。`--label`等は元の集計と同じ値を指定します。集計から行う場合は`--plots-only`を外します。各解析のオプションはそれぞれのREADMEを参照してください。

`--result-dir`と`--log-dir`を指定するときは、保存したい解析名のdirectoryまで含めたパスを渡します。

## 移動済みの既存出力

従来`analysis/result/channel_response/`直下にあった応答比較の出力は、その下の`adc_energy/`に移動しました。従来のジョブログ・task manifest・投入記録も`job/adc_energy/`に移動しました。ファイル内容は変更していません。

過去のログ、設定記録、ROOT内の入力由来情報には、当時の実行パスがそのまま残っています。再実行には上記の現行コマンドを使ってください。

旧パス用の互換symlinkと転送Makefileは、投入済みジョブの終了確認後に削除しました。新規実行・投入は各解析のsubdirectoryから行います。
