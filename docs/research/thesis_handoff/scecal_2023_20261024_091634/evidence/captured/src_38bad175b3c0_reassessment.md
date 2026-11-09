ScECAL beam test data / simulation 不一致の再検討（2026-09-29）

添付の前回回答に加え、保存済み比較ROOT、100 GeVのcalibration ROOTの標本、対応するMC raw・truth、および現在配置されている校正・digitizationソースを確認した。既存の解析ソース・校正定数・ROOTは変更していない。

今回の結論は、**MIP校正の代替定数に具体的な不具合があり、それとは別に、縦方向の信号分布、HG/LGの整合性、周辺hitの再現にも問題が見える**というもの。中央channelのhit数の一致は応答の一致を保証しない。現時点で、単独の原因が全差分を説明したとは言えない。

今回作成した[比較図](/home/murata_t/scecal_reassessment_20260929/reassessment_summary.png)、[PDF](/home/murata_t/scecal_reassessment_20260929/reassessment_summary.pdf)、[保存fit値の一覧](/home/murata_t/scecal_reassessment_20260929/saved_fit_values.csv)、[100 GeVの層別平均](/home/murata_t/scecal_reassessment_20260929/saved_layer_means_100GeV.csv)を参照。

1. 新たに確認したMIP校正の問題

現在の `Calibration.cxx` の `MIPFactorExtract_simulation()` は、10 µm channelで `ChiSquare/NDF <= 2` を満たすMIP値を平均するが、MIP値が正かどうかを検査していない。

現在の設定で参照される `all_auto_muon_v4_trackfit.root` を読むと、10 µm側でこの条件を通る3381件のうち490件が非正値であり、最小値は -10 ADCだった。これらを含める平均は **81.576853 ADC**、正値に限定した同条件の平均は **97.098354 ADC** になる。

fit品質条件で除外されたchannelには前者の平均が代入される。したがって、平均値への異常値混入が、正常な信号を測定しているchannelのエネルギー換算にも波及する。

同じADC信号を97.10で割る代わりに81.58で割れば、換算値は約19.0%大きくなる。ただし、97.10をただちに正しい代替定数とみなすこともできない。まず失敗値を除外した上で、channel別の校正の妥当性を確認する必要がある。

さらにdigitizationの `Extract::MIPPeakValueExtract()` では、負のMIP値を除外する一方、上記のχ²/NDFカットを使用していない。このため**同じ定数ファイルを渡しても、信号生成に使うMIP値と、ADCからエネルギーへ戻すMIP値が一致しない**。

100 GeVの標本でHG/LGを調べられた中央10 µmの185 channelのうち、160 channelが校正側のχ²/NDF条件に該当して代替値を使う対象だった。これは全検出器に対する割合ではない。

ソースの問題だけで終わらせず、保存済みMCのrawとcalibrationの同一イベント・同一CellIDを対応させ、温度が20℃であることを確認して

\[
M_{\rm inferred}=0.305\,(ADC_{\rm HG}-pedestal)/E_{\rm HG}
\]

を計算したところ、複数channelで **81.57685 ADC** が再現した。例えばCellID 420008は定数ファイルでは87.5975 ADCだが、保存ROOTからの逆算値は81.57685 ADCである。この問題は、現在のソースに存在するだけでなく、少なくとも調べた保存済みMCに現れている。

ただし、この不具合の修正だけでdata/MC比が1になるとは断定できない。dataにも同じ代替値が使われていれば、共通の影響があるためである。各channelについて、dataの実応答、MC生成時の定数、両者の校正時の定数を対応させる必要がある。

根拠：[Calibration.cxx:369](/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src/Calibration.cxx:369)、[Extract.cxx:186](/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/src/Extract.cxx:186)、[監査集計](/home/murata_t/scecal_reassessment_20260929/mip_audit_summary.json)、[保存MCから逆算した定数](/home/murata_t/scecal_reassessment_20260929/mip_inferred_from_saved_mc.csv)。

2. 中央イベントを選び直してもenergy差は残る

既存の `analysis/result/cog_channel_scan/` に前回提案のCoG×channel比較が保存されていた。`beforeEventCut` の保存Gaussian中心値を読み直した結果は以下。比はdata/MC、channel範囲はx・yの半幅22.5 mm、CoGはchannelカット前の全有効hitから計算したlayer 9・10の両方に対する条件。

| Beam energy | CoGなし：energy比 | CoG 5 mm：energy比 | CoG 5 mm：Nhit比 |
|---|---:|---:|---:|
| 30 GeV | 1.0378 | 1.0999 | 0.9841 |
| 50 GeV | 1.0676 | 1.1282 | 0.9890 |
| 100 GeV | 1.0839 | 1.1476 | 0.9964 |
| 150 GeV | 1.0971 | 1.1006 | 1.0021 |
| 250 GeV | 1.1216 | 1.1374 | 1.0067 |

100 GeV、CoG 5 mm、中央channelでNhitを200–204に限定した2次元ヒストグラムの射影では、平均energyはdata **4087.45 MeV**、MC **3581.30 MeV**、比は **1.1413**。したがって、別々にfitしたNhit分布とenergy分布を比べたことだけで生じた見かけの差ではない。

ただし同じNhitであることは、同じshower形状・真のエネルギー沈着量・入射位置を保証しない。Nhit固定だけで校正原因に一意に絞ることはできない。

また、CoGを5 mmに固定するとenergy比は単調増加ではなくなる。前回の「高エネルギーほど差が大きいので大信号非線形性」という推論は、選択条件・run依存性を除いてから評価する必要がある。

根拠：[条件付き平均の数値](/home/murata_t/scecal_reassessment_20260929/conditional_energy_100GeV.csv)、[解析でのCoG選択](/home/murata_t/ScECAL_BeamTest/analysis/cog_channel_scan/SingleEnergyAnalysis.cc:451)。fit中心値の統計誤差はCSVに保存したが、fit範囲・校正などの系統誤差を含まない。

3. 縦方向の分布差を、前回より優先して調べるべき

100 GeV、CoG 5 mm、channel半幅22.5 mmでの層別平均energyは、data/MC比がlayer 5で1.359、9で1.224、10で1.167。一方、layer 20で0.847、24で0.723、28で0.748になる。layer番号は0始まり。

層別平均energyから作る重心 \(\sum_l l\langle E_l\rangle/\sum_l\langle E_l\rangle\) はdata **12.178層**、MC **12.813層**。再構成された信号はdataの方が約0.64層前に寄っている。channel半幅68.2 mmでもdata12.528層、MC13.221層で、同方向の差がある。

これは真のshower開始位置の測定結果ではない。層ごとの校正差でもこの形を作れる。候補は、MIP代替値・HG/LG・SiPM応答の層依存、上流物質、absorber厚さや材料・gap、入射角、縦方向のleakageなどである。まず校正の層依存を調べ、その後にtruthのshower開始位置・縦方向分布を比較するのが妥当。

一様なenergy倍率だけでこの前後の差を消すことはできない。layer 14のような局所的な落ち込みは、滑らかなshower形状の差とは分けてchannel・layer校正を確認したい。

中央のlayer 9で8 hitとなる割合はdata95.4%、MC95.5%。hit数の頭打ちは実際に存在するが、それはenergy差を見えにくくする条件であり、energy差そのものを発生させる原因ではない。

4. HG/LGのforward変換と逆変換が整合していない

現在のMCはpedestalを除くと、おおむね \(LG=s\,HG\) として低gain信号を作る。一方、校正側はfitされた切片bを使って \(HG_{\rm reco}=(LG-b)/s\) と戻す。MC生成側に同じbがなければ、\(-b/s\) の余分なoffsetが生じる。

参照されるinter-calibrationファイルでは、10 µmの有効slope範囲にあるchannelの切片中央値は約-4.06 LG ADC、HG換算offset中央値は約+145 ADCだった。

100 GeVのdata Run105・106・107から各先頭2000イベント、MCファイル0から1000イベントを確認した。bad mask適用後の中央10 µm channel、保存HG energyが2–5 MeVの範囲で、\(E_{LG}-E_{HG}\) のhit単位中央値はdataで約-0.01〜+0.04 MeV、MCで **+0.845 MeV**。同じchannelに限定してchannelごとの中央値を比較しても、MC-dataの中央値は約+0.85〜+0.88 MeVだった。

この確認はHG/LGの重なり領域の診断であり、2–5 MeVでLG枝が最終energyとして採用されているという意味ではない。観測されたLG側の正offsetは、LGが選ばれるとMC energyを高くする方向に働く。したがって**dataの方が高い今回の差を単独では説明せず、別の差を打ち消している可能性がある**。修正するとdata/MC差が大きくなる可能性もあるが、変換の整合性は先に確保すべき。

根拠：[MCのHG/LG生成](/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/src/MCDigi.cxx:264)、[校正の切片](/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src/Calibration.cxx:453)、[標本の集計](/home/murata_t/scecal_reassessment_20260929/calibration_sample_100GeV.csv)。

5. SiPM飽和・crosstalkは、MIPでの正規化をそろえて検証する

MCでは実測MIP ADC / SiPM gainで光量/MIPを決め、その後10 µmでは10%、15 µmでは15%のcrosstalkを追加している。実測MIP電荷にすでに含まれるcrosstalkまで一次光電子として扱えば二重計上になる。現在のソースには補正らしき行があるがコメントアウトされている。

この効果は低信号でMC電荷・threshold通過数を増やす方向で、周辺MC hitの過剰に寄与し得る。中央のdata energy過剰の単独説明にはならない。MIP値がgainと光量の両方に使われるため、単純にcrosstalkだけOFFにして総energyを合わせる検証は不十分。まず低信号・MIPの応答を同じ基準に合わせ、その状態で大信号の曲がり方を比較する。

対応するMC truthとcalibrationの先頭500イベントを比較すると、10 µmのtruth 2–5 MeVで再構成energy / truthのhit単位平均は約1.300だった。これはthreshold通過後の標本で、全channel共通の校正倍率ではない。さらに中央channelのraw/calib/truthを用いた先頭50イベントの確認では、truth 1–3 MeVに対するrawの低信号応答は、ファイル中MIP値を使った線形期待の約1.105倍（中央値）。ファイルMIP / 逆算された校正MIPは中央値約1.202だった。MIP変換の不整合とcrosstalk追加が重なり得ることを示している。

一方、`PixelCorrection()` は1回のchannel応答につき各pixelを一度しか発火させず、時間的回復・再発火を明示的に扱わない。実機で光パルス中の再発火が効けば、MCが実機より強く飽和して、大信号でMC energyが低くなる候補になる。ただし実機の光パルス幅・回復時間を確認する必要がある。メーカーも飽和応答がpixel数だけでなく回復時間と光パルス幅に依存すると説明している。[浜松ホトニクスの技術資料](https://www.hamamatsu.com/us/en/resources/interactive-tools/mppc-sipm-linearity.html)

根拠：[MIPから光量への変換とcrosstalk](/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/src/MCDigi.cxx:216)、[pixel発火モデル](/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/src/MCDigi.cxx:140)、[MC truth比較](/home/murata_t/scecal_reassessment_20260929/mc_truth_calibration_closure_sample.csv)。

6. 周辺hit、beam条件、温度などの残る候補

| 候補 | 今回の根拠・限界 | 切り分け |
|---|---|---|
| threshold・効率・周辺noise | 100 GeV、CoG 5 mmでも68.2 mmのNhit比は0.894、中央は0.996。周辺に差が残る。ただしMCのshowerが広い可能性もある | channel別の小信号効率と、shower軸からの距離別Nhitを比較。飽和・gainを調整する前に低信号のMIP基準を一致させる |
| beam位置・幅・角度、strip内の集光位置 | CoG条件でenergy比が変わる。CoG 5 mm通過率は100 GeVでdata8.4%・MC11.7%、150 GeVでdata31.0%・MC12.0%と変わる | 外部trackerがあれば入射位置・角度をそろえる。なければrun別にCoGと幅を比較。通過率だけからbeam幅を逆算しない |
| stripの長手方向の光量応答 | channel中心座標のカットは真の入射位置のカットではない。MCの4.2%一様乱数は位置に相関した応答の再現とは異なる | strip長手位置・層の向き別にMIPとshower信号を比較 |
| run依存のgain・温度・bias | 100 GeVの標本ではdataの温度中央値22.86–23.09℃、MC20℃。同条件の中央energy中央値は3runで約4092–4141 MeV、MC3584 MeVであり、標本内では単一runだけの異常ではない | 全runで温度とMIP基準温度・biasの対応を確認。100 GeV以外もrun別に調べる |
| 温度再構成の未初期化変数 | 現在の`tempReconstruction()`は`localTemperature`を初期化せず加算する。コードとして不具合 | 修正対象だが、標本の保存温度は非物理的には見えない。今回の主要因と断定しない。符号・係数・基準温度も別途確認 |
| contamination・pileup・fit方法 | beforeEventCutであり、CoGなしの30 GeVなどはfit比と全平均比が大きく異なる | 同じイベント選択、同じbin幅で、fit中心・中央値・条件付き平均を併記。E–Nhitカットで差を人工的に消さない |
| ADC非線形・HG切替・clipping・timing | 現在の校正はraw HG 2600で切替、MCのHG上限判定は別のchannel別値 | raw HG–LG相関、切替前後、LG大信号の傾きを確認。HGの切替とSiPM pixel飽和を分離 |
| dead channel、DAQ・memory cell、event alignment | 共通bad maskだけではrun固有の欠損や重複を保証できない | CellIDの物理channel化、重複数、layer別の欠損・時刻・memory cell依存を点検 |
| geometry・上流物質・EM model・quenching | 縦方向と周辺分布の差の候補だが、今回Geant4設定を新規に検証したわけではない | 校正の整合を確保後、absorber・支持材・gap・上流物質、truthの開始位置・漏れを優先。physics listやproduction cutはその後 |

温度の根拠：[未初期化変数](/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src/EBUdecode.cxx:70)、[補正式](/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/src/Calibration.cxx:639)。

7. 次の検証はこの順序がよい

1. **使用定数の監査とMIP代替値の修正検証。** data校正・MC生成・MC校正について、各CellIDの採用MIP値、元fit値、採用/不採用理由、代替元、定数ファイルを保存する。負値・非有限値を代替平均に混ぜない。dataとMCの生成・校正での実効定数を比較する。修正後は同じrawに対する差分を比較する。
2. **HG/LGの変換を閉じる。** 同じslope・intercept・pedestal基準で、線形入力→digitization→calibrationが期待応答に戻るかをchannel別に検証する。HG/LG重なり領域で見えたMCのoffsetを解消する。
3. **同じtruthイベントで応答段階を分ける。** 線形応答、光電子統計、crosstalk、SiPM飽和、ADC、thresholdを順に加え、各段階のE・Nhit・layer分布を保存する。低信号/MIPの正規化をそろえてから飽和の曲率を比べる。
4. **run・位置をそろえた層別比較。** 100 GeVを起点に、CoG 5 mmと22.5/68.2 mmの両範囲で、E_l、N_l、中央/広範囲のenergy割合を比較する。校正の層差が消えてもdataが前方寄りなら、上流物質・geometry・shower開始位置の優先度を上げる。
5. **周辺の低信号効率を独立に検証。** 同一イベントで中心と周辺を分解し、しきい値より十分上のhitとthreshold付近を分ける。記録されなかったdata hitは後からthresholdを下げても戻せないため、後段の一律カットだけで効率曲線を決めない。

12 MeVで打ち切られる既存1-hitヒストグラムは、大信号一致の判定には使えない。100 GeV、CoG 5 mm、中央10 µmではoverflowはdata40.2%、MC40.9%だった。割合が近くても、overflow内の電荷の大きさは不明である。全範囲の枝から、hit数分布とエネルギー寄与を別々に作る必要がある。

今回の比較は保存結果の再読と限定したイベント標本の監査であり、全rawの再校正・全MCの再digitizationは行っていない。現在の全ソースが保存ROOTの生成時と同一だったかも完全には追跡していない。ただしMIP代替値81.57685 ADCは保存MCから直接再現でき、HG/LGのずれと中央のenergy差は保存ROOT自体で確認できた。

再現用のスクリプトは `/home/murata_t/scecal_reassessment_20260929/read_saved_results.py`、`/home/murata_t/scecal_reassessment_20260929/probe_calibration.py`、`/home/murata_t/scecal_reassessment_20260929/audit_mip_constants.py`、`/home/murata_t/scecal_reassessment_20260929/plot_summary.py`。実行環境はROOT対応の `python3.12`。結果はスクリプトのある `/home/murata_t/scecal_reassessment_20260929/` に出力する。

資料の参照先: [索引](README.md)、[関連する図・表・実行記録](/home/murata_t/scecal_reassessment_20260929/)。
