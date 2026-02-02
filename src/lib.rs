use rtsyn_plugin::{PluginApi, PluginString};
use serde_json::Value;
use std::ffi::c_void;

const INPUTS: &[&str] = &["i_syn"];
const OUTPUTS: &[&str] = &["x", "y", "z"];

#[derive(Debug)]
struct HindmarshRosev2Rust {
    x: f64,
    y: f64,
    z: f64,
    input_syn: f64,
    e: f64,
    mu: f64,
    s: f64,
    vh: f64,
    dt: f64,
    burst_duration: f64,
    s_points: usize,
    period_seconds: f64,
    cfg_x: f64,
    cfg_y: f64,
    cfg_z: f64,
}

impl HindmarshRosev2Rust {
    fn new() -> Self {
        let x = -0.9013747551021072;
        let y = -3.15948829665501;
        let z = 3.247826955037619;
        Self {
            x,
            y,
            z,
            input_syn: 0.0,
            e: 3.25,
            mu: 0.006,
            s: 4.0,
            vh: 1.0,
            dt: 0.0015,
            burst_duration: 1.0,
            s_points: 1,
            period_seconds: 0.001,
            cfg_x: x,
            cfg_y: y,
            cfg_z: z,
        }
    }

    fn update_burst_settings(&mut self) {
        if self.period_seconds <= 0.0 {
            self.s_points = 1;
            return;
        }
        let freq = 1.0 / self.period_seconds;
        if self.burst_duration > 0.0 {
            let pts_burst = set_pts_burst(self.burst_duration, freq, &mut self.dt);
            let mut s_points = (pts_burst / (self.burst_duration * freq)).round() as usize;
            if s_points == 0 {
                s_points = 1;
            }
            self.s_points = s_points;
        } else {
            let steps = ((self.period_seconds / self.dt).round() as usize).max(1);
            self.s_points = steps;
        }
        if self.s_points == 0 {
            self.s_points = 1;
        }
    }

    fn set_config(&mut self, config: &Value) {
        let get = |key: &str, default: f64| -> f64 {
            config.get(key).and_then(|v| v.as_f64()).unwrap_or(default)
        };
        let x = get("x", self.x);
        let y = get("y", self.y);
        let z = get("z", self.z);
        if (x, y, z) != (self.cfg_x, self.cfg_y, self.cfg_z) {
            self.cfg_x = x;
            self.cfg_y = y;
            self.cfg_z = z;
            self.x = x;
            self.y = y;
            self.z = z;
        }
        self.e = get("e", self.e);
        self.mu = get("mu", self.mu);
        self.s = get("s", self.s);
        self.vh = get("vh", self.vh);
        self.dt = get("time_increment", self.dt).max(0.0);
        self.burst_duration = get("burst_duration", self.burst_duration);
        self.period_seconds = get("period_seconds", self.period_seconds);
        self.update_burst_settings();
    }

    fn process(&mut self) {
        let dt = self.dt;
        let steps = self.s_points.min(10_000).max(1);
        for _ in 0..steps {
            let mut vars = [self.x, self.y, self.z];
            let mut k = [[0.0f64; 3]; 6];
            let mut aux = [0.0f64; 3];

            let f = |vars: [f64; 3], params: &Self| -> [f64; 3] {
                let x = vars[0];
                let y = vars[1];
                let z = vars[2];
                let v =
                    y + 3.0 * (x * x) - (x * x * x) - params.vh * z + params.e - params.input_syn;
                let ydot = 1.0 - 5.0 * (x * x) - y;
                let zdot = params.mu * (-params.vh * z + params.s * (x + 1.6));
                [v, ydot, zdot]
            };

            let r0 = f(vars, self);
            for j in 0..3 {
                k[0][j] = dt * r0[j];
                aux[j] = vars[j] + k[0][j] * 0.2;
            }

            let r1 = f(aux, self);
            for j in 0..3 {
                k[1][j] = dt * r1[j];
                aux[j] = vars[j] + k[0][j] * 0.075 + k[1][j] * 0.225;
            }

            let r2 = f(aux, self);
            for j in 0..3 {
                k[2][j] = dt * r2[j];
                aux[j] = vars[j] + k[0][j] * 0.3 - k[1][j] * 0.9 + k[2][j] * 1.2;
            }

            let r3 = f(aux, self);
            for j in 0..3 {
                k[3][j] = dt * r3[j];
                aux[j] =
                    vars[j] + k[0][j] * 0.075 + k[1][j] * 0.675 - k[2][j] * 0.6 + k[3][j] * 0.75;
            }

            let r4 = f(aux, self);
            for j in 0..3 {
                k[4][j] = dt * r4[j];
                aux[j] = vars[j] + k[0][j] * 0.660493827160493 + k[1][j] * 2.5
                    - k[2][j] * 5.185185185185185
                    + k[3][j] * 3.888888888888889
                    - k[4][j] * 0.864197530864197;
            }

            let r5 = f(aux, self);
            for j in 0..3 {
                k[5][j] = dt * r5[j];
            }

            for j in 0..3 {
                vars[j] += k[0][j] * 0.098765432098765
                    + k[2][j] * 0.396825396825396
                    + k[3][j] * 0.231481481481481
                    + k[4][j] * 0.308641975308641
                    - k[5][j] * 0.035714285714285;
            }

            self.x = vars[0];
            self.y = vars[1];
            self.z = vars[2];
        }
    }
}

fn set_pts_burst(burst_duration: f64, freq: f64, dt: &mut f64) -> f64 {
    let dts: [f64; 144] = [
        0.000500, 0.000600, 0.000700, 0.000800, 0.000900, 0.001000, 0.001100, 0.001200, 0.001300,
        0.001400, 0.001500, 0.001600, 0.001800, 0.002000, 0.002200, 0.002500, 0.002800, 0.002900,
        0.003000, 0.003100, 0.003200, 0.003300, 0.003400, 0.003500, 0.003600, 0.003700, 0.003800,
        0.003900, 0.004000, 0.004100, 0.004200, 0.004300, 0.004400, 0.004500, 0.004600, 0.004700,
        0.004800, 0.004900, 0.005000, 0.005100, 0.005200, 0.005400, 0.005600, 0.005800, 0.006000,
        0.006200, 0.006400, 0.006600, 0.006800, 0.007000, 0.007200, 0.007400, 0.007700, 0.008000,
        0.008300, 0.008600, 0.008900, 0.009200, 0.009600, 0.010000, 0.010400, 0.010900, 0.011400,
        0.011900, 0.012500, 0.013100, 0.013800, 0.014600, 0.015400, 0.016300, 0.017300, 0.018500,
        0.019900, 0.021500, 0.023300, 0.025500, 0.028100, 0.028400, 0.028700, 0.029000, 0.029400,
        0.029800, 0.030200, 0.030600, 0.031000, 0.031400, 0.031800, 0.032200, 0.032600, 0.033000,
        0.033400, 0.033900, 0.034400, 0.034900, 0.035400, 0.035900, 0.036400, 0.036900, 0.037400,
        0.038000, 0.038600, 0.039200, 0.039800, 0.040400, 0.041000, 0.041700, 0.042400, 0.043100,
        0.043800, 0.044500, 0.045300, 0.046100, 0.046900, 0.047700, 0.048600, 0.049500, 0.050400,
        0.051400, 0.052400, 0.053400, 0.054500, 0.055600, 0.056800, 0.058000, 0.059300, 0.060600,
        0.062000, 0.063400, 0.064900, 0.066500, 0.068200, 0.069900, 0.071700, 0.073600, 0.075600,
        0.077700, 0.079900, 0.082300, 0.084800, 0.087500, 0.090300, 0.093300, 0.096500, 0.100000,
    ];
    let pts: [f64; 144] = [
        577638.000000,
        481366.000000,
        412599.000000,
        357615.500000,
        317880.000000,
        286092.500000,
        259143.333333,
        237548.000000,
        218869.500000,
        203236.000000,
        189687.000000,
        177634.000000,
        157897.000000,
        142001.833333,
        129024.142857,
        113496.125000,
        101304.555556,
        97811.222222,
        94527.400000,
        91478.200000,
        88619.400000,
        85916.636364,
        83389.636364,
        81007.090909,
        78743.583333,
        76615.416667,
        74599.250000,
        72676.000000,
        70859.076923,
        69130.846154,
        67476.642857,
        65907.357143,
        64402.666667,
        62971.466667,
        61602.533333,
        60286.187500,
        59030.250000,
        57825.562500,
        56664.411765,
        55553.294118,
        54485.000000,
        52463.222222,
        50586.263158,
        48841.842105,
        47211.050000,
        45685.666667,
        44255.818182,
        42914.772727,
        41650.739130,
        40459.083333,
        39335.208333,
        38270.680000,
        36778.346154,
        35398.000000,
        34117.571429,
        32926.517241,
        31815.833333,
        30777.612903,
        29493.939394,
        28313.588235,
        27223.638889,
        25974.405405,
        24834.410256,
        23790.268293,
        22647.767442,
        21609.977778,
        20513.166667,
        19388.627451,
        18381.132075,
        17365.719298,
        16361.600000,
        15299.937500,
        14223.202899,
        13164.400000,
        12147.123457,
        11098.876404,
        10071.693878,
        9965.282828,
        9861.100000,
        9759.059406,
        9626.242718,
        9497.009615,
        9371.179245,
        9248.672897,
        9129.293578,
        9012.981818,
        8899.594595,
        8789.000000,
        8681.149123,
        8575.896552,
        8473.170940,
        8348.176471,
        8226.809917,
        8108.934426,
        7994.379032,
        7883.015873,
        7774.710938,
        7669.348837,
        7566.801527,
        7447.308271,
        7331.525926,
        7219.291971,
        7110.435714,
        7004.816901,
        6902.298611,
        6786.417808,
        6674.355705,
        6565.940397,
        6460.987013,
        6359.339744,
        6247.012579,
        6138.592593,
        6033.866667,
        5932.660714,
        5822.777778,
        5716.896552,
        5614.796610,
        5505.541436,
        5400.467391,
        5299.324468,
        5192.348958,
        5089.615385,
        4982.075000,
        4878.985294,
        4772.014354,
        4669.633803,
        4564.178899,
        4463.381166,
        4360.214912,
        4255.294872,
        4149.216667,
        4048.296748,
        3946.654762,
        3844.764479,
        3743.041353,
        3641.872263,
        3541.583630,
        3438.300000,
        3336.926421,
        3233.951299,
        3133.666667,
        3032.899696,
        2932.320588,
        2829.684659,
    ];
    select_dt_neuron_model(&dts, &pts, burst_duration * freq, dt)
}

fn select_dt_neuron_model(dts: &[f64], pts: &[f64], pts_live: f64, dt: &mut f64) -> f64 {
    let mut aux = pts_live;
    let mut factor = 1.0;
    let mut pts_burst = -1.0;
    let mut selected_dt = -1.0;
    let mut flag = false;

    while aux < pts[0] {
        aux = pts_live * factor;
        factor += 1.0;
        for i in (0..pts.len()).rev() {
            if pts[i] > aux {
                selected_dt = dts[i];
                pts_burst = pts[i];
                let fract = (pts_burst / pts_live).fract();
                let intp = (pts_burst / pts_live).trunc();
                if fract <= 0.1 * intp {
                    flag = true;
                }
                break;
            }
        }
        if flag {
            break;
        }
    }

    if !flag {
        for i in (0..pts.len()).rev() {
            if pts[i] > aux {
                selected_dt = dts[i];
                pts_burst = pts[i];
                break;
            }
        }
    }

    if selected_dt > 0.0 {
        *dt = selected_dt;
    }
    pts_burst
}

extern "C" fn create(_id: u64) -> *mut c_void {
    let instance = Box::new(HindmarshRosev2Rust::new());
    Box::into_raw(instance) as *mut c_void
}

extern "C" fn destroy(handle: *mut c_void) {
    if handle.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(handle as *mut HindmarshRosev2Rust));
    }
}

extern "C" fn meta_json(_handle: *mut c_void) -> PluginString {
    let value = serde_json::json!({
        "name": "Hindmarsh Rose Dyn Rust",
        "kind": "hindmarsh_rose_dyn_rs"
    });
    PluginString::from_string(value.to_string())
}

extern "C" fn inputs_json(_handle: *mut c_void) -> PluginString {
    PluginString::from_string(serde_json::to_string(INPUTS).unwrap_or_default())
}

extern "C" fn outputs_json(_handle: *mut c_void) -> PluginString {
    PluginString::from_string(serde_json::to_string(OUTPUTS).unwrap_or_default())
}

extern "C" fn set_config_json(handle: *mut c_void, data: *const u8, len: usize) {
    if handle.is_null() || data.is_null() || len == 0 {
        return;
    }
    let slice = unsafe { std::slice::from_raw_parts(data, len) };
    if let Ok(json) = serde_json::from_slice::<Value>(slice) {
        let instance = unsafe { &mut *(handle as *mut HindmarshRosev2Rust) };
        instance.set_config(&json);
    }
}

extern "C" fn set_input(handle: *mut c_void, name: *const u8, len: usize, value: f64) {
    if handle.is_null() || name.is_null() || len == 0 {
        return;
    }
    let slice = unsafe { std::slice::from_raw_parts(name, len) };
    if let Ok(name) = std::str::from_utf8(slice) {
        if name == "i_syn" {
            let instance = unsafe { &mut *(handle as *mut HindmarshRosev2Rust) };
            instance.input_syn = value;
        }
    }
}

extern "C" fn process(handle: *mut c_void, _tick: u64) {
    if handle.is_null() {
        return;
    }
    let instance = unsafe { &mut *(handle as *mut HindmarshRosev2Rust) };
    instance.process();
}

extern "C" fn get_output(handle: *mut c_void, name: *const u8, len: usize) -> f64 {
    if handle.is_null() || name.is_null() || len == 0 {
        return 0.0;
    }
    let slice = unsafe { std::slice::from_raw_parts(name, len) };
    if let Ok(name) = std::str::from_utf8(slice) {
        let instance = unsafe { &mut *(handle as *mut HindmarshRosev2Rust) };
        return match name {
            "x" => instance.x,
            "y" => instance.y,
            "z" => instance.z,
            _ => 0.0,
        };
    }
    0.0
}

#[no_mangle]
pub extern "C" fn rtsyn_plugin_api() -> *const PluginApi {
    static API: PluginApi = PluginApi {
        create,
        destroy,
        meta_json,
        inputs_json,
        outputs_json,
        set_config_json,
        set_input,
        process,
        get_output,
    };
    &API as *const PluginApi
}
