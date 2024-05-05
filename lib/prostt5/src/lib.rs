#[cfg(feature = "mkl")]
extern crate intel_mkl_src;

#[cfg(feature = "accelerate")]
extern crate accelerate_src;

use anyhow::Result;
use candle_core::utils::{cuda_is_available, metal_is_available};
use candle_core::{DType, Device, IndexOp, Tensor};
use candle_nn::{ops::softmax, VarBuilder};
use candle_transformers::models::t5::{self, T5EncoderModel};
use serde_json;
use std::collections::HashMap;
// use std::fs::File;
// use std::io::Write;
use std::path::PathBuf;

pub mod cnn;
use crate::cnn::CNN;

pub fn device(cpu: bool) -> Result<Device> {
    if cpu {
        Ok(Device::Cpu)
    } else if cuda_is_available() {
        Ok(Device::new_cuda(0)?)
    } else if metal_is_available() {
        Ok(Device::new_metal(0)?)
    } else {
        Ok(Device::Cpu)
    }
}

const DTYPE: DType = DType::F16;
const BIT_FACTOR: f32 = 8.0;
const SCORE_BIAS: f32 = 0.0;
const PROFILE_AA_SIZE: i32 = 20;

fn compute_log_pssm(profile: Vec<f32>, query_length: usize) -> Result<Vec<i8>> {
    let pback: Vec<f32> = vec![
        0.0489372, 0.0306991, 0.101049, 0.0329671, 0.0276149, 0.0416262, 0.0452521, 0.030876,
        0.0297251, 0.0607036, 0.0150238, 0.0215826, 0.0783843, 0.0512926, 0.0264886, 0.0610702,
        0.0201311, 0.215998, 0.0310265, 0.0295417, 0.00001,
    ];
    //let pback: Vec<f32> = vec![0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.052, 0.00001];
    let mut pssm: Vec<i8> = vec![0; query_length * PROFILE_AA_SIZE as usize];
    println!("{query_length}");
    println!("{}", profile.len());
    for pos in 0..query_length {
        for aa in 0..PROFILE_AA_SIZE {
            let aa_prob: f32 = profile[pos * PROFILE_AA_SIZE as usize + aa as usize];
            let idx: usize = pos * PROFILE_AA_SIZE as usize + aa as usize;
            let log_prob: f32 = (aa_prob / pback[aa as usize]).log2();
            let mut pssmval: f32 = log_prob * BIT_FACTOR + BIT_FACTOR * SCORE_BIAS;
            if pssmval < 0.0 {
                pssmval -= 0.5;
            } else {
                pssmval += 0.5;
            }
            if pssmval > 127.0 {
                pssmval = 127.0;
            }
            if pssmval < -128.0 {
                pssmval = -128.0;
            }
            pssm[idx] = pssmval as i8;
        }
    }
    Ok(pssm)
}

fn to_buffer(consensus: &[u8], center_sequence_len: usize, pssm: Vec<i8>) -> Result<Vec<u8>> {
    let pssm_u8 = unsafe { std::slice::from_raw_parts(pssm.as_ptr() as *const u8, pssm.len()) };
    // for prob in pssm_u8 {
    //     print!("{prob}\n");
    // }
    let mut result: Vec<u8> = Vec::new();
    for pos in 0..center_sequence_len - 1 {
        for aa in 0..PROFILE_AA_SIZE {
            result.push(pssm_u8[pos * PROFILE_AA_SIZE as usize + aa as usize]);
        }
        result.push(char_to_number(consensus[pos]));
        result.push(char_to_number(consensus[pos]));
        result.push(0 as u8);
        result.push(0 as u8);
        result.push(0 as u8);
    }
    result.push(0 as u8);
    Ok(result)
}

fn number_to_char(n: u32) -> u8 {
    match n {
        0 => b'A',
        1 => b'C',
        2 => b'D',
        3 => b'E',
        4 => b'F',
        5 => b'G',
        6 => b'H',
        7 => b'I',
        8 => b'K',
        9 => b'L',
        10 => b'M',
        11 => b'N',
        12 => b'P',
        13 => b'Q',
        14 => b'R',
        15 => b'S',
        16 => b'T',
        17 => b'V',
        18 => b'W',
        19 => b'Y',
        _ =>  b'X', // Default case for numbers not in the list
    }
}

fn char_to_number(n: u8) -> u8 {
    match n {
        b'A' => 0,
        b'C' => 1,
        b'D' => 2,
        b'E' => 3,
        b'F' => 4,
        b'G' => 5,
        b'H' => 6,
        b'I' => 7,
        b'K' => 8,
        b'L' => 9,
        b'M' => 10,
        b'N' => 11,
        b'P' => 12,
        b'Q' => 13,
        b'R' => 14,
        b'S' => 15,
        b'T' => 16,
        b'V' => 17,
        b'W' => 18,
        b'Y' => 19,
        _ => 20, // Default case for numbers not in the list
    }
}

pub struct T5ModelBuilder {
    pub profile: bool,
    pub device: Device,
    pub config: t5::Config,
    pub weights_filename: Vec<PathBuf>,
    pub cnn_filename: Vec<PathBuf>,
    pub tokens_map: HashMap<String, usize>,
}

impl T5ModelBuilder {
    pub fn load(base_path: &str, profile: bool, cpu: bool, cache: bool) -> Result<Self> {
        let device = device(cpu)?;

        let base = PathBuf::from(base_path);
        let config_filename = base.join("config.json");

        let weights_filename = vec![base.join("model.safetensors")];
        let path = if profile {
            base.join("new_cnn.safetensors")
        } else {
            base.join("cnn.safetensors")
        };
        let cnn_filename = vec![path];
        let config = std::fs::read_to_string(config_filename)?;
        let mut config: t5::Config = serde_json::from_str(&config)?;
        config.use_cache = cache;

        let tokens_filename = base.join("tokens.json");
        let tokens_config = std::fs::read_to_string(tokens_filename)?;
        let tokens_map: HashMap<String, usize> = serde_json::from_str(&tokens_config)?;

        Ok(Self {
            profile,
            device,
            config,
            weights_filename,
            cnn_filename,
            tokens_map,
        })
    }

    pub fn build_encoder(&self) -> Result<T5EncoderModel> {
        let vb = unsafe {
            VarBuilder::from_mmaped_safetensors(&self.weights_filename, DTYPE, &self.device)?
        };
        Ok(T5EncoderModel::load(vb, &self.config)?)
    }

    pub fn build_cnn(&self) -> Result<CNN> {
        let vb = unsafe {
            VarBuilder::from_mmaped_safetensors(&self.cnn_filename, DTYPE, &self.device)?
        };
        Ok(CNN::load(vb, self.profile)?)
    }
}

pub struct ProstT5 {
    pub builder: T5ModelBuilder,
    pub encoder: T5EncoderModel,
    pub cnn: CNN,
}

impl ProstT5 {
    pub fn load(base_path: String, profile: bool, cpu: bool, cache: bool) -> Result<ProstT5> {
        let builder = T5ModelBuilder::load(
            &base_path,
            profile,
            cpu,
            cache
        )?;
        let encoder = builder.build_encoder()?;
        let cnn = builder.build_cnn()?;

        Ok(Self {
            builder,
            encoder,
            cnn
        })
    }

    pub fn predict(
        &mut self,
        sequence: String,
        // buffer_file: &mut File,
        // index_file: &mut File,
        // buffer_file_seqs: &mut File,
        // index_file_seqs: &mut File,
        // embeds_file: &mut File,
        // index: i32,
        // curr_len: i32,
        // seq_len: i32,
    ) -> Result<(Vec<u8>, i32)> {
        let copy = &sequence;
        // Replace each character in the string
        let replaced_values: Vec<Option<&usize>> = copy
            .chars()
            .map(|c| self.builder.tokens_map.get(&c.to_string()))
            .collect();

        let unknown_value: usize = 2; // Default value for None

        let ts: Vec<&usize> = replaced_values
            .iter()
            .map(|option| option.unwrap_or(&unknown_value))
            .collect();
        let mut tokens: Vec<&usize> = vec![self.builder.tokens_map.get("<AA2fold>").unwrap()];
        tokens.extend(ts.iter().clone());
        tokens.push(self.builder.tokens_map.get("</s>").unwrap());
        let tokens: Vec<i64> = tokens.iter().map(|&num| *num as i64).collect();

        let input_token_ids = Tensor::new(&tokens[..], &self.builder.device)?
            .unsqueeze(0)?
            .to_dtype(DType::U8)?;

        let ys = self.encoder.forward(&input_token_ids)?;

        let ys = ys.i((.., 1..ys.dims3()?.1 - 1))?;
        let ys = ys.pad_with_zeros(1, 0, 1)?;
        // let _ = writeln!(embeds_file, "{}", ys);

        let output = &self.cnn.forward(&ys)?;

        let vals = output.argmax_keepdim(1)?;
        let vals = vals.flatten(0, 2)?;
        let vals = vals.to_vec1::<u32>()?;
        let structure: Vec<u8> = vals.iter().map(|&n| number_to_char(n)).collect();

        if self.builder.profile {
            let probs: Tensor = softmax(output, 1)?;
            let lines: Vec<Vec<f32>> = probs.to_dtype(DType::F32)?.to_vec3()?[0].to_vec();

            let mut profile: Vec<f32> = Vec::new();
            for i in 0..lines[0].len() {
                for j in 0..20 {
                    profile.push(lines[j][i]);
                }
            }

            let consensus = structure.clone();
            let mut pssm = compute_log_pssm(profile, structure.len())?;
            for pos in 0..consensus.len() {
                if consensus[pos] == b'X' {
                    for aa in 0..PROFILE_AA_SIZE {
                        pssm[pos * PROFILE_AA_SIZE as usize + aa as usize] = -1;
                    }
                }
            }
            let result = to_buffer(&consensus, structure.len(), pssm)?;
            let length = result.len() as i32;
            // buffer_file.write(&result)?;

            // let index_content = format!("{}\t{}\t{}\n", index, curr_len, result.len());
            // index_file.write(index_content.as_bytes())?;

            // write!(buffer_file_seqs, "{}\n\0", &prompt)?;

            // let index_content = format!("{}\t{}\t{}\n", index, seq_len, &prompt.len() + 2);
            // index_file_seqs.write(index_content.as_bytes())?;

            Ok((result, length))
        } else {
            Ok((structure, 0))
        }
    }
}
