# Bindgen Cuda

[![Latest version](https://img.shields.io/crates/v/bindgen_cuda.svg)](https://crates.io/crates/bindgen_cuda)
[![Documentation](https://docs.rs/bindgen_cuda/badge.svg)](https://docs.rs/bindgen_cuda)
![License](https://img.shields.io/crates/l/bindgen_cuda.svg)

Similar crate than [bindgen](https://github.com/rust-lang/rust-bindgen) in philosophy.
It will help create automatic bindgen to cuda kernels source files and make them easier to use
directly from Rust.


## PTX inclusion
Let's say you have a file

`src/cuda.cu`
```cuda
__global__ void cuda_hello(){
    printf("Hello World from GPU!\n");
}
```

You can add `bindgen_cuda` as a build dependency:

```bash
cargo add --build bindgen_cuda
```

And then create this `build.rs`

```no_run
fn main() {
    let builder = bindgen_cuda::Builder::default();
    let bindings = builder.build_ptx().unwrap();
    bindings.write("src/lib.rs");
}
```

This will create a src file containing the following code:

```ignore
pub const CUDA: &str = include_str!(concat!(env!("OUT_DIR"), "/cuda.ptx"));
```

You can then use the PTX directly in your rust code with a library like [cudarc](https://github.com/coreylowman/cudarc/).

## Raw cuda calls
Alternatively you can build a static library that you can link against in build.rs in order to call cuda directly with the c code.

`src/cuda.cu`

```cuda
__global__ void cuda_hello(){
    printf("Hello World from GPU!\n");
}

int run() {
    cuda_hello<<<1,1>>>(); 
    return 0;
}
```


Then write the `build.rs`:

```no_run
fn main() {
    let builder = bindgen_cuda::Builder::default();
    builder.build_lib("libcuda.a");
    println!("cargo:rustc-link-lib=cuda");
}
```

Which you can then interface through FFI in `src/lib.rs`:


```no_run
extern "C" {
    fn cuda_hello();
}
fn main(){
    unsafe{ cuda_hello();}
}
```

