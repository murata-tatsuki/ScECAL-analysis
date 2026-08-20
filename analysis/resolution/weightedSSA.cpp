#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// 最終的な出力を格納する構造体
struct SSAResult {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> layer;
    std::vector<double> energy;
};

// 内部処理用のストリップ構造体
struct StripHit {
    double x;
    double y;
    double energy;
};

// 1次元の重なり（交差長さ）を高速に計算するインライン関数
inline double GetOverlap1D(double a1, double a2, double b1, double b2) {
    double overlap = std::min(a2, b2) - std::max(a1, b1);
    return (overlap > 0.0) ? overlap : 0.0;
}

/**
 * @brief 重なり面積(Overlap Area)を用いた重み付きSSAを実行する関数（ノイズ除去機能付き）
 * * @param in_x        入力: 各ヒットのX中心座標 (mm)
 * @param in_y        入力: 各ヒットのY中心座標 (mm)
 * @param in_layer    入力: 各ヒットのレイヤー番号 (0.0 ~ 29.0)
 * @param in_energy   入力: 各ヒットのエネルギー
 * @return SSAResult  分割後の仮想セル情報を含んだ構造体
 */
SSAResult PerformStripSplitSSA(const std::vector<double>& in_x,
                               const std::vector<double>& in_y,
                               const std::vector<double>& in_layer,
                               const std::vector<double>& in_energy) 
{
    SSAResult result;
    
    // データの整合性チェック
    if (in_x.empty() || in_x.size() != in_y.size() || 
        in_x.size() != in_layer.size() || in_x.size() != in_energy.size()) {
        return result; 
    }

    // 【最適化1】1イベントのヒットをレイヤーごと(0~29)にバケットソート形式で整理
    // これにより、上下レイヤーの探索が O(N) から O(1) レベルに高速化されます
    std::vector<std::vector<StripHit>> layer_hits(30);
    for (size_t i = 0; i < in_x.size(); ++i) {
        int lay = static_cast<int>(in_layer[i] + 0.5); 
        if (lay >= 0 && lay < 30) {
            layer_hits[lay].push_back({in_x[i], in_y[i], in_energy[i]});
        }
    }

    // 【最適化2】42個配置されている「細かなピッチベース」の中心座標を予めルックアップテーブル化
    std::vector<double> C42(42);
    for (int i = 0; i < 42; ++i) {
        C42[i] = -108.65 + i * 5.3;
    }

    // 物理形状パラメータ (ミリ単位)
    const double L_long = 45.0;            // ストリップの長軸長さ
    const double L_short = 5.0;            // ストリップの短軸幅
    const double H_long = L_long / 2.0;    // 半長 (22.5 mm)
    const double H_short = L_short / 2.0;  // 半幅 (2.5 mm)

    // レイヤーごとのループ処理
    for (int lay = 0; lay < 30; ++lay) {
        if (layer_hits[lay].empty()) continue;

        bool is_even = (lay % 2 == 0);

        // 重み付けに使用する上下のレイヤーIDを特定
        std::vector<int> adj_layers;
        if (lay - 1 >= 0) adj_layers.push_back(lay - 1);
        if (lay + 1 < 30) adj_layers.push_back(lay + 1);

        // 各ストリップヒットの処理
        for (const auto& hit : layer_hits[lay]) {
            if (hit.energy <= 0) continue;

            if (is_even) {
                // ==========================================
                // 偶数レイヤー: X方向に長く(45mm)、Y方向に短い(5mm)
                // ➡ X方向を 42個のグリッドへ分割する
                // ==========================================
                std::vector<double> weights(42, 0.0);
                std::vector<double> base_areas(42, 0.0);
                double total_weight = 0.0;

                for (int k = 0; k < 42; ++k) {
                    double X_v = C42[k];
                    // 仮想セルが、この長細いストリップの物理境界内に収まっているか（重なり）
                    double v_x1 = std::max(X_v - H_short, hit.x - H_long);
                    double v_x2 = std::min(X_v + H_short, hit.x + H_long);
                    double b_area = (v_x2 > v_x1) ? (v_x2 - v_x1) * L_short : 0.0;
                    
                    if (b_area <= 0.0) continue; // 幾何学的に全く重なっていない
                    base_areas[k] = b_area;

                    // 上下レイヤー(奇数レイヤー)のエネルギーとの「2次元重なり面積」を算出
                    double w_k = 0.0;
                    for (int alay : adj_layers) {
                        for (const auto& ahit : layer_hits[alay]) {
                            double ol_x = GetOverlap1D(v_x1, v_x2, ahit.x - H_short, ahit.x + H_short);
                            double ol_y = GetOverlap1D(hit.y - H_short, hit.y + H_short, ahit.y - H_long, ahit.y + H_long);
                            w_k += ahit.energy * (ol_x * ol_y);
                        }
                    }
                    weights[k] = w_k;
                    total_weight += w_k;
                }

                // 【条件変更】上下レイヤーにヒットがない（total_weight == 0）場合は、ノイズとみなして除外
                if (total_weight <= 0.0) {
                    continue; 
                }

                // 重みに基づいてエネルギーを分配し、構造体に格納
                for (int k = 0; k < 42; ++k) {
                    if (base_areas[k] <= 0.0 || weights[k] <= 0.0) continue;
                    
                    double final_w = weights[k] / total_weight;
                    result.x.push_back(C42[k]);
                    result.y.push_back(hit.y); // Y方向は元から42グリッドの細かさ
                    result.layer.push_back(static_cast<double>(lay));
                    result.energy.push_back(hit.energy * final_w);
                }
            } 
            else {
                // ==========================================
                // 奇数レイヤー: X方向に短く(5mm)、Y方向に長い(45mm)
                // ➡ Y方向を 42個のグリッドへ分割する
                // ==========================================
                std::vector<double> weights(42, 0.0);
                std::vector<double> base_areas(42, 0.0);
                double total_weight = 0.0;

                for (int m = 0; m < 42; ++m) {
                    double Y_v = C42[m];
                    // 仮想セルのY物理境界と、この縦長ストリップの境界の重なり
                    double v_y1 = std::max(Y_v - H_short, hit.y - H_long);
                    double v_y2 = std::min(Y_v + H_short, hit.y + H_long);
                    double b_area = (v_y2 > v_y1) ? L_short * (v_y2 - v_y1) : 0.0;

                    if (b_area <= 0.0) continue;
                    base_areas[m] = b_area;

                    // 上下レイヤー(偶数レイヤー)のエネルギーとの「2次元重なり面積」を算出
                    double w_m = 0.0;
                    for (int alay : adj_layers) {
                        for (const auto& ahit : layer_hits[alay]) {
                            double ol_x = GetOverlap1D(hit.x - H_short, hit.x + H_short, ahit.x - H_long, ahit.x + H_long);
                            double ol_y = GetOverlap1D(v_y1, v_y2, ahit.y - H_short, ahit.y + H_short);
                            w_m += ahit.energy * (ol_x * ol_y);
                        }
                    }
                    weights[m] = w_m;
                    total_weight += w_m;
                }

                // 【条件変更】上下レイヤーにヒットがない（total_weight == 0）場合は、ノイズとみなして除外
                if (total_weight <= 0.0) {
                    continue; 
                }

                // エネルギーの分配・構造体への格納
                for (int m = 0; m < 42; ++m) {
                    if (base_areas[m] <= 0.0 || weights[m] <= 0.0) continue;

                    double final_w = weights[m] / total_weight;
                    result.x.push_back(hit.x); // X方向は元から42グリッドの細かさ
                    result.y.push_back(C42[m]);
                    result.layer.push_back(static_cast<double>(lay));
                    result.energy.push_back(hit.energy * final_w);
                }
            }
        }
    }

    return result;
}