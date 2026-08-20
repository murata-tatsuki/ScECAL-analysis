#!/bin/bash

cd /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/e-/sps
# cd /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/sps

# すべての *GeV ディレクトリを順に処理
for energy_dir in *GeV; do
    # 念のためディレクトリであることを確認
    if [ -d "$energy_dir" ]; then
        
        # 中にある threshold* ディレクトリを順に処理
        for thresh_path in "$energy_dir"/threshold*; do
            if [ -d "$thresh_path" ]; then
                # 末尾のディレクトリ名だけを取得 (例: threshold10)
                thresh_name=$(basename "$thresh_path")
                
                # 移動先の親ディレクトリ(threshold*)を作成
                mkdir -p "$thresh_name"
                
                # ディレクトリごと新しい階層に移動
                # 例: 0.5GeV/threshold10 -> threshold10/0.5GeV
                mv "$thresh_path" "$thresh_name/$energy_dir"
            fi
        done
        
        # threshold* をすべて移動し終えて空になった *GeV ディレクトリを削除
        # （他にファイルが残っている場合は安全のため削除されません）
        rmdir "$energy_dir" 2>/dev/null
    fi
done

echo "ディレクトリ構造の変換が完了しました！"