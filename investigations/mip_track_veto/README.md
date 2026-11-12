# 貫通する MIP 型 track の veto の簡易調査

5、40、100 GeV の data・単電子 MC に、疎で細い貫通 track の形状条件を適用する。
正式な解析マクロや較正を変更せず、結果は `investigations/result/mip_track_veto/` に保存する。

## 入力と判定

前回の `multiple_shower_veto` investigation の各 sample/energy 2万事象をそのまま使用する。
sampling manifest、event ID、信号和、CoG、data の HG–LG tail 除外、現行 MC の再監査済み tail flag を引き継ぐ。
前回の複数 shower veto 自体は適用しない。入力309ファイルの size/mtime が前回記録と一致することを確認する。

nominal は `analysis/shower_leak/TrackSelection.hh` の v1 の幾何条件を使用し、次のすべてを要求する。

- `0 < good-hit 信号和 <= 0.03 GeV`（ビーム energy ではなく校正済み信号）
- 採用 hit 数60以下、正の層信号がある active layer 数15以上
- 最初の5層、最後の5層でそれぞれ active layer 数2以上
- hit 数 / active layer 数が2以下
- 奇数層 X、偶数層 Y にそれぞれ4層以上の信号がある
- 各投影の幾何 RMS が5.3 mm以下

層内の CoG・幅は energy 重み付き、層間の直線 fit は各層に等しい重みの OLS。
幾何 RMS は `sqrt(mean((layer CoG - fitted line)^2 + within-layer width^2))` とする。
層位置は現在の EBUdecode の値を使用し、保存 `Hit_Z` との一致も確認する。
hit 数には採用されたゼロ energy hit と重複 memory hit を個別に数える。新しい hit threshold は加えない。

信号和・mask は、前回と同じ**元の full CellID による mask lookup**を維持する。
既存 `shower_leak/Aggregate.cc` の physical CellID lookup へ変更しないため、既存解析全体の完全再現ではない。
同じ基準で before/after を比較する調査であり、physical-ID mask での信号和も診断値として保存する。

CoG は layer 9・10 の各 X/Y の矩形条件で200/20/10/5 mm。
200 mm の図表示は既存規則に合わせ `nocut` とするが、両層の正の信号条件は残る。
そのため `E<=0.03 GeV` かつ200 mm条件を通る事象だけを現 `Calib_Hit` から再読すれば、4条件すべての veto を判定できる。
再読した全候補について ID と前回信号和の一致を要求する。高信号事象の track を探索した調査ではない。

必要層数20への強化、幾何 RMS 10.6 mmへの緩和も事前固定した参考条件として保存する。
この tag は μ/非相互作用πと整合する形状の選別であり、粒子同定や真の混入率の測定ではない。

## 幅と解釈

前回と同じ処方で、全体の RMS/mean、中央68%半幅/median、簡易 Gaussian core σ/μ を保存する。
histogram の bin と before で定めた core fit 範囲を after にも固定する。
core の範囲が30 MeVより上なら、今回の veto はその範囲内の事象を除かないので、core幅が不変になる。
その場合でも低信号成分・平均・全体の幅は変わり得るため、RMS/meanと分けて報告する。
Gaussian fit は正式な resolution マクロと同じ処方ではない。

nominal の全体幅・68%幅の差には、同じ event を再抽出する400回の paired bootstrap 区間を付ける。
区間は採用 sample 内の事象統計だけで、run/spill相関、sampling、較正・モデルの系統誤差は含まない。

## 再現と出力

プロジェクトのルートから実行する。前回結果が必要。既存の完了済み結果を再生成する場合は別の `--out` を指定する。

```bash
python3 investigations/mip_track_veto/check.py --energy 5
python3 investigations/mip_track_veto/check.py --energy 40
python3 investigations/mip_track_veto/check.py --energy 100
python3 investigations/mip_track_veto/check.py --combine
```

- `comparison/summary20mm.csv`：代表値の一覧
- `comparison/{metrics,tag_rates,paired_changes}.csv`：全energy・CoG・判定条件の数値
- `comparison/rms_summary.png`：全体 RMS/mean の比較
- `{data,mc}/diagnostics/<energy>GeV/`：event特徴、energy分布、候補track図
- `<sample>/<condition>/<CoG>mm/<energy>GeV/`：条件別の数値・histogram
- `inputs/`：元のsampling manifest、信号和・入力検証、実行条件・使用ソース

幅の計算には前回の `multiple_shower_veto/check.py` の関数を読み取り利用し、使用版を各 `inputs` に保存する。

## 相互作用pionを狙う追加調査

`pion_shower.py` は、上記の貫通MIP vetoに「前方のMIP型trackから途中でshowerへ移る事象」のtagを追加する。
出力は同じ調査結果の `pion_shower/` 以下。電子data・単電子MCは前回と同じ各2万事象、同じtail除外・good hit信号和・CoG条件を使う。
旧貫通veto、追加shower tag、それらのunionを別々に保存する。正式な解析コードは変更しない。

追加tagの事前固定条件は以下のとおり。層番号は0始まり。

- shower開始候補 k=8–25を探索。looseは6–25、tightは10–25として開始下限だけを変える。
- layer 0..k−1の75%以上がactiveで、奇数X・偶数Yそれぞれ3層以上に信号がある。
- 前方active層の信号中央値が0.3–2.5 MIP、4 MIPを超える層は最大1層、正信号の固有strip数/active層数が2.5以下。
- 前方のX・Yの層内幅込み直線RMSが各5.3 mm以下。fitの層間重みは等しい。
- kとk+1の連続2層で、各層信号が6 MIP以上、正信号の固有stripが各4本以上。
- k..k+3の平均層信号が、前方active層の平均信号の5倍以上。

1 MIP=0.305 MeV。同じphysical stripのmemory重複は、新tagのstrip数では1本として数える。
前方の信号和・重心・幅には前回と同じ採用hitを使う。**追加tagに総event信号の上限・energy窓は設けない**。
これらの閾値をresolutionの改善率で調整しない。
早期に相互作用するπや、最後の層だけでshowerを始めるπは取り逃す。遅れてshowerを作る電子やμの放射もtagされ得る。
単電子MCの除去率と幅変化を同時に確認し、πの確定PIDや混入率の測定とは区別する。

π設定で取得された `calib_original/pi-` からも各energy6千事象を全fileに比例配分して形状参考に使う。
PID情報・対応tail判定・現行電子較正との生成条件一致が確認できないため、このsampleには分解能を出さず、候補率を純πの除去効率と呼ばない。
旧π較正の信号scale差も参考候補率に影響し得る。

```bash
python3 investigations/mip_track_veto/pion_shower.py --energy 5
python3 investigations/mip_track_veto/pion_shower.py --energy 40
python3 investigations/mip_track_veto/pion_shower.py --energy 100
python3 investigations/mip_track_veto/pion_shower.py --combine
```

`pion_shower/comparison/summary20mm.csv` が代表値。全CoG条件の数値・図、追加tag例、入力manifest・検証・ソースを同じ結果配下に保持する。
`pion_shower/pion_reference/data/` はπ設定runの形状参考で、電子data/MCの幅比較とは分ける。
`--summarize-only` は保存済み層別momentから同じ固定条件で再集計する場合に用いる。

### 早期showerを要求するπ抑制用の選別

前方MIP prefixを厳密に要求する条件は、π設定runでも相互作用候補の取り逃しが多かった。
そのため `--selection early` では、前方でshowerが発達する電子型事象を残す固定条件を別途比較する。
結果は `pion_shower/early_shower/` に保存し、上記prefix条件の結果は保持する。

- 層信号6 MIP以上、正信号の固有physical strip4本以上を満たす層をdenseとする。
- 連続2層がdenseとなる最初の先頭層kをshower開始の指標とする。
- nominalは **k<=6を保持**。層番号は0始まりなので、7番目の層までに始まり、次の層もdenseであることを要求する。
- denseな連続2層がどこにもない事象は除外する。
- 既存の貫通MIP vetoもORで維持する。looseはk<=8、tightはk<=4を保持する。
- 総event信号のenergy窓や、前方MIP prefixの条件は新しい早期shower選別には用いない。

層別momentは上の調査で検証した同一事象から読み取り再利用し、入力元のSHAを記録する。
π設定参考sampleの候補率と単電子MCの誤除去率を見て、この条件をresolutionの改善量を見る前に固定した。
5 GeVでは遅いshower発達の単電子も除くので、dataの幅改善をπ除去だけの効果と解釈しない。
早期に相互作用するπは残る。この選別はπ抑制用であり、πの確定PIDではない。

```bash
python3 investigations/mip_track_veto/pion_shower.py --selection early --energy 5
python3 investigations/mip_track_veto/pion_shower.py --selection early --energy 40
python3 investigations/mip_track_veto/pion_shower.py --selection early --energy 100
python3 investigations/mip_track_veto/pion_shower.py --selection early --combine
```

`early_shower/comparison/summary20mm.csv` はdata/MCの候補率とbefore/after幅、既存MIP vetoだけの場合の幅も含む。
`early_shower/comparison/pion_reference_rates.csv` の値は較正・PID・tail選択が一致しないπ設定runでの参考値に限る。
