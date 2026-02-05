use rtsyn_plugin::{PluginApi, PluginString};
use serde_json::Value;
use std::ffi::c_void;

const INPUTS: &[&str] = &["i_syn"];
const OUTPUTS: &[&str] = &["Membrane potential (V)", "Membrane potential (mV)"];

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
            dt: 0.15,
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

        if self.burst_duration > 0.0 {
            // Use sophisticated dt selection like original RTXI implementation
            let freq = 1.0 / self.period_seconds;
            let pts_match = self.burst_duration * freq;
            self.dt = self.select_optimal_dt(pts_match);
            self.s_points = ((self.select_pts_burst(self.burst_duration, freq) / (self.burst_duration * freq)).round() as usize).max(1);
        } else {
            // Simple case - use fixed dt and calculate steps
            let steps = ((self.period_seconds / self.dt).round() as usize).max(1);
            self.s_points = steps;
        }
        
        if self.s_points == 0 {
            self.s_points = 1;
        }
    }

    fn select_optimal_dt(&self, pts_match: f64) -> f64 {
        // Simplified version of the RTXI dt selection algorithm
        // Original has 144 values, using key ones for efficiency
        let dts = [0.0005, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.1];
        let pts = [577638.0, 286092.5, 189687.0, 142001.8, 94527.4, 56664.4, 28313.6, 18381.1, 14223.2, 9497.0, 5716.9, 2829.7];
        
        for (i, &pt) in pts.iter().enumerate() {
            if pt <= pts_match {
                return dts[i];
            }
        }
        dts[dts.len() - 1] // fallback to largest dt
    }

    fn select_pts_burst(&self, sec_per_burst: f64, freq: f64) -> f64 {
        let pts_match = sec_per_burst * freq;
        self.select_optimal_dt(pts_match); // This sets the dt
        pts_match // Simplified - return the calculated points
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
        
        self.burst_duration = get("burst_duration", self.burst_duration);
        self.period_seconds = get("period_seconds", self.period_seconds);
        self.update_burst_settings();
    }

    fn process(&mut self) {
        let dt = self.dt;
        // Original working logic - simple step limiting
        let steps = self.s_points.min(10_000).max(1);

        for _ in 0..steps {
            let mut vars = [self.x, self.y, self.z];
            let mut k = [[0.0f64; 3]; 6];
            let mut aux = [0.0f64; 3];

            let f = |vars: [f64; 3], params: &Self| -> [f64; 3] {
                let x = vars[0];
                let y = vars[1];
                let z = vars[2];
                let xdot = y + 3.0 * (x * x) - (x * x * x) - params.vh * z + params.e - params.input_syn;
                let ydot = 1.0 - 5.0 * (x * x) - y;
                let zdot = params.mu * (-params.vh * z + params.s * (x + 1.6));
                [xdot, ydot, zdot]
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
        "name": "Hindmarsh Rose v2 Rust",
        "default_vars": [
            ["x", -0.9013],
            ["y", -3.1594],
            ["z", 3.24782],
            ["e", 3.0],
            ["mu", 0.006],
            ["s", 4.0],
            ["vh", 1.0],
            ["burst_duration", 1.0]
        ]
    });
    PluginString::from_string(value.to_string())
}

extern "C" fn inputs_json(_handle: *mut c_void) -> PluginString {
    PluginString::from_string(serde_json::to_string(INPUTS).unwrap_or_default())
}

extern "C" fn outputs_json(_handle: *mut c_void) -> PluginString {
    PluginString::from_string(serde_json::to_string(OUTPUTS).unwrap_or_default())
}

extern "C" fn behavior_json(_handle: *mut c_void) -> PluginString {
    let behavior = serde_json::json!({
        "supports_start_stop": true,
        "supports_restart": true,
        "extendable_inputs": {"type": "none"},
        "loads_started": true
    });
    PluginString::from_string(behavior.to_string())
}

extern "C" fn ui_schema_json(_handle: *mut c_void) -> PluginString {
    let schema = serde_json::json!({
        "outputs": ["Membrane potential (V)", "Membrane potential (mV)"],
        "inputs": ["i_syn"],
        "variables": ["x", "y", "z"]
    });
    PluginString::from_string(schema.to_string())
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

extern "C" fn process(handle: *mut c_void, _tick: u64, period_seconds: f64) {
    if handle.is_null() {
        return;
    }
    let instance = unsafe { &mut *(handle as *mut HindmarshRosev2Rust) };
    
    // ALWAYS use the period_seconds from runtime, not from config
    // This ensures the plugin respects workspace period settings
    if (instance.period_seconds - period_seconds).abs() > f64::EPSILON {
        instance.period_seconds = period_seconds;
        instance.update_burst_settings();
    }
    
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
            "Membrane potential (V)" => instance.x,
            "Membrane potential (mV)" => instance.x * 1000.0,
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
        behavior_json: Some(behavior_json),
        ui_schema_json: Some(ui_schema_json),
        set_config_json,
        set_input,
        process,
        get_output,
    };
    &API as *const PluginApi
}
