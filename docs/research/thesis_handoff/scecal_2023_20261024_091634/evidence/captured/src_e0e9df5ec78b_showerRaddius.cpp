#include <vector>
#include <cmath>
#include <algorithm>

// レイヤーごとの解析結果を格納する構造体
struct LayerShower {
    double cog_x = 0.0;       // そのレイヤーの重心 X
    double cog_y = 0.0;       // そのレイヤーの重心 Y
    double rms_radius = 0.0;  // そのレイヤーの RMS半径
    double r90_radius = 0.0;  // そのレイヤーの R90半径
    double total_energy = 0.0;// そのレイヤーの 総堆積エネルギー
};

// 全体の解析結果を格納する構造体（レイヤーごとの結果も内包）
struct ShowerContainer {
    double cog_x = 0.0;
    double cog_y = 0.0;
    double rms_radius = 0.0;
    double r90_radius = 0.0;
    std::vector<LayerShower> layers; // サイズ30の配列として使用
};

// R90のソート用補助構造体
struct CellDistance {
    double radius;
    double energy;
    bool operator<(const CellDistance& other) const {
        return radius < other.radius;
    }
};

/**
 * @brief 全体およびレイヤーごとのシャワー重心と半径(RMS, R90)を一度に計算する関数
 * @param ssa_x        SSA通過後のx座標
 * @param ssa_y        SSA通過後のy座標
 * @param ssa_layer    SSA通過後のレイヤー番号 (追加)
 * @param ssa_energy   SSA通過後のenergy
 * @param fraction     包含半径の比率 (デフォルト0.90)
 * @return ShowerContainer 全体とレイヤーごとの計算結果を一括で返す
 */
ShowerContainer CalculateShowerRadius(const std::vector<double>& ssa_x,
                                      const std::vector<double>& ssa_y,
                                      const std::vector<double>& ssa_layer,
                                      const std::vector<double>& ssa_energy,
                                      double fraction = 0.90) 
{
    ShowerContainer shower;
    shower.layers.resize(30); // 0~29レイヤー分の枠を確保

    size_t n_cells = ssa_energy.size();
    if (n_cells == 0 || ssa_x.size() != n_cells || ssa_y.size() != n_cells || ssa_layer.size() != n_cells) {
        return shower;
    }

    // =============================================================
    // Pass 1: 全体と各レイヤーの総エネルギー ＆ 重心(CoG) を同時に計算
    // =============================================================
    double global_total_energy = 0.0;
    double global_sum_xe = 0.0;
    double global_sum_ye = 0.0;

    std::vector<double> layer_sum_xe(30, 0.0);
    std::vector<double> layer_sum_ye(30, 0.0);

    for (size_t i = 0; i < n_cells; ++i) {
        double e = ssa_energy[i];
        if (e <= 0.0) continue;
        
        int lay = static_cast<int>(ssa_layer[i] + 0.5);
        if (lay < 0 || lay >= 30) continue;

        // 全体の加算
        global_total_energy += e;
        global_sum_xe += ssa_x[i] * e;
        global_sum_ye += ssa_y[i] * e;

        // レイヤーごとの加算
        shower.layers[lay].total_energy += e;
        layer_sum_xe[lay] += ssa_x[i] * e;
        layer_sum_ye[lay] += ssa_y[i] * e;
    }

    if (global_total_energy <= 0.0) return shower;

    // 重心の確定
    shower.cog_x = global_sum_xe / global_total_energy;
    shower.cog_y = global_sum_ye / global_total_energy;

    for (int lay = 0; lay < 30; ++lay) {
        if (shower.layers[lay].total_energy > 0.0) {
            shower.layers[lay].cog_x = layer_sum_xe[lay] / shower.layers[lay].total_energy;
            shower.layers[lay].cog_y = layer_sum_ye[lay] / shower.layers[lay].total_energy;
        }
    }

    // =============================================================
    // Pass 2: RMS半径の計算 ＆ R90計算用の距離リスト作成（全体＋各レイヤー）
    // =============================================================
    double global_sum_r2_e = 0.0;
    std::vector<CellDistance> global_cell_list;
    global_cell_list.reserve(n_cells);

    std::vector<double> layer_sum_r2_e(30, 0.0);
    std::vector<std::vector<CellDistance>> layer_cell_lists(30);

    for (size_t i = 0; i < n_cells; ++i) {
        double e = ssa_energy[i];
        if (e <= 0.0) continue;

        int lay = static_cast<int>(ssa_layer[i] + 0.5);
        if (lay < 0 || lay >= 30) continue;

        // 全体の重心からの距離
        double dx_g = ssa_x[i] - shower.cog_x;
        double dy_g = ssa_y[i] - shower.cog_y;
        double r_g = std::sqrt(dx_g * dx_g + dy_g * dy_g);
        
        global_sum_r2_e += (r_g * r_g) * e;
        global_cell_list.push_back({r_g, e});

        // そのレイヤーの重心からの距離
        double dx_l = ssa_x[i] - shower.layers[lay].cog_x;
        double dy_l = ssa_y[i] - shower.layers[lay].cog_y;
        double r_l = std::sqrt(dx_l * dx_l + dy_l * dy_l);

        layer_sum_r2_e[lay] += (r_l * r_l) * e;
        layer_cell_lists[lay].push_back({r_l, e});
    }

    // RMS半径の確定
    shower.rms_radius = std::sqrt(global_sum_r2_e / global_total_energy);
    for (int lay = 0; lay < 30; ++lay) {
        if (shower.layers[lay].total_energy > 0.0) {
            shower.layers[lay].rms_radius = std::sqrt(layer_sum_r2_e[lay] / shower.layers[lay].total_energy);
        }
    }

    // =============================================================
    // Pass 3: 距離でソートし、包含半径 (R90) を特定する
    // =============================================================
    
    // 全体の R90 計算
    std::sort(global_cell_list.begin(), global_cell_list.end());
    double target_energy = global_total_energy * fraction;
    double accumulated_energy = 0.0;
    for (const auto& cell : global_cell_list) {
        accumulated_energy += cell.energy;
        if (accumulated_energy >= target_energy) {
            shower.r90_radius = cell.radius;
            break;
        }
    }

    // レイヤーごとの R90 計算
    for (int lay = 0; lay < 30; ++lay) {
        if (shower.layers[lay].total_energy <= 0.0) continue;

        std::sort(layer_cell_lists[lay].begin(), layer_cell_lists[lay].end());
        double target_layer_energy = shower.layers[lay].total_energy * fraction;
        double accumulated_layer_energy = 0.0;
        
        for (const auto& cell : layer_cell_lists[lay]) {
            accumulated_layer_energy += cell.energy;
            if (accumulated_layer_energy >= target_layer_energy) {
                shower.layers[lay].r90_radius = cell.radius;
                break;
            }
        }
    }

    return shower;
}