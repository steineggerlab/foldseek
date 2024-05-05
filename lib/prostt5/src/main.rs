use anyhow::Result;
use pico_args::Arguments;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use std::path::Path;
use prostt5::ProstT5;
pub mod cnn;

#[cfg(feature = "tracing")]
use tracing_chrome::ChromeLayerBuilder;
#[cfg(feature = "tracing")]
use tracing_subscriber::prelude::*;

struct Args {
    cpu: bool,
    disable_cache: bool,
    prompt: Option<String>,
    generate_profile: bool,
    output: Option<String>,
}

fn process_fasta(
    input_path: &Path,
    output_path: &Path,
    prost: &mut ProstT5
) -> io::Result<()> {
    let file = File::open(input_path)?;
    let reader = io::BufReader::new(file);

    let mut output_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(output_path)?;
    let mut time_file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("./times.fasta")?;
    // let mut buffer_file = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile_ss")?;
    // let mut index_file = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile_ss.index")?;
    // let mut dbtype_file = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile_ss.dbtype")?;
    // let mut buffer_file_seqs = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile")?;
    // let mut index_file_seqs = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile.index")?;
    // let mut dbtype_file_seqs = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile.dbtype")?;
    // let mut buffer_file_h = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile_h")?;
    // let mut index_file_h = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile_h.index")?;
    // let mut dbtype_file_h = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile_h.dbtype")?;
    // let mut lookup_file = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./profile.lookup")?;
    // let mut embeds_file = OpenOptions::new()
    //     .create(true)
    //     .write(true)
    //     .truncate(true)
    //     .open("./embeds.txt")?;

    // let u8_vec_ss: Vec<u8> = vec![2, 0, 0, 0];
    // let u8_vec_aa: Vec<u8> = vec![0, 0, 0, 0];
    // let u8_vec_h: Vec<u8> = vec![12, 0, 0, 0];
    // dbtype_file.write(&u8_vec_ss)?;
    // dbtype_file_seqs.write(&u8_vec_aa)?;
    // dbtype_file_h.write(&u8_vec_h)?;
    let mut s = String::new();
    // let mut index = 0;
    // let mut curr_len = 0;
    // let mut seq_len = 0;
    // let mut header_len = 0;

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            // let l = line.clone();
            // write!(buffer_file_h, "{}\n\0", &l[1..])?;
            // let index_content = format!("{}\t{}\t{}\n", index, header_len, l.len() + 1);
            // index_file_h.write(index_content.as_bytes())?;
            // header_len += (l.len() + 1) as i32;
            // let index_content = format!("{}\t{}\t{}\n", index, &l[1..], 0);
            // lookup_file.write(index_content.as_bytes())?;
            if !s.is_empty() {
                let start_time = std::time::Instant::now();
                let prompt = s.clone();
                let prediction = prost.predict(
                    prompt,
                    // &mut buffer_file,
                    // &mut index_file,
                    // &mut buffer_file_seqs,
                    // &mut index_file_seqs,
                    // &mut embeds_file,
                    // index,
                    // curr_len,
                    // seq_len,
                );
                let prediction = prediction.unwrap().0;

                // curr_len += prediction.1;
                // seq_len += (s.len() + 2) as i32;
                println!("Took {:?}", start_time.elapsed());
                unsafe {
                    let _ = writeln!(output_file, "{}", std::str::from_utf8_unchecked(&prediction));
                }
                let _ = writeln!(
                    time_file,
                    "{} : {}",
                    s.len(),
                    start_time.elapsed().as_secs() as f64
                        + start_time.elapsed().subsec_nanos() as f64 * 1e-9
                );
                s.clear();
            }
            writeln!(output_file, "{}", line)?;
        } else {
            // index += 1;
            s.push_str(&line);
        }
    }

    // Write the last sequence if not empty
    if !s.is_empty() {
        let prompt = s.clone();
        let prediction = prost.predict(
            prompt,
            // &mut buffer_file,
            // &mut index_file,
            // &mut buffer_file_seqs,
            // &mut index_file_seqs,
            // &mut embeds_file,
            // index,
            // curr_len,
            // seq_len,
        );
        let res = prediction.unwrap().0;
        unsafe {
            let _ = writeln!(output_file, "{}", std::str::from_utf8_unchecked(&res));
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let mut args = Arguments::from_env();

    // Convert the argument parsing to manually handle each option
    let cpu = args.contains("--cpu");
    let disable_cache = args
        .opt_value_from_str("--disable-cache")
        .unwrap_or(Some(false))
        .unwrap_or(false);
    let prompt = args.opt_value_from_str("--prompt").unwrap_or(None);
    let generate_profile = args
        .opt_value_from_str("--generate-profile")
        .unwrap_or(Some(false))
        .unwrap_or(false);
    let output = args.opt_value_from_str("--output").unwrap_or(None);

    // Construct Args
    let args = Args {
        cpu,
        disable_cache,
        prompt,
        generate_profile,
        output,
    };

    #[cfg(feature = "tracing")]
    let _guard = {
        let (chrome_layer, guard) = ChromeLayerBuilder::new().build();
        tracing_subscriber::registry().with(chrome_layer).init();
        Some(guard)
    };

    let mut prost = ProstT5::load("model/".to_string(), args.generate_profile, args.cpu, !args.disable_cache)?;

    match args.prompt {
        Some(prompt) => {
            let output = args.output.unwrap();
            let start = std::time::Instant::now();
            let _ = process_fasta(
                Path::new(&prompt),
                Path::new(&output),
                &mut prost,
            );
            println!("Took {:?}", start.elapsed());
        }
        None => {}
    }
    Ok(())
}
