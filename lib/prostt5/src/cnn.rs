use anyhow::Result;
use candle_core::{Tensor, D};
use candle_nn::{ops::log_softmax, Conv2d, Conv2dConfig, Module, VarBuilder};

pub fn conv2d_non_square(
    in_channels: usize,
    out_channels: usize,
    kernel_size1: usize,
    kernel_size2: usize,
    cfg: Conv2dConfig,
    vb: VarBuilder,
) -> Result<Conv2d> {
    let init_ws = candle_nn::init::DEFAULT_KAIMING_NORMAL;
    let ws = vb.get_with_hints(
        (
            out_channels,
            in_channels / cfg.groups,
            kernel_size1,
            kernel_size2,
        ),
        "weight",
        init_ws,
    )?;
    let bound = 1. / (in_channels as f64).sqrt();
    let init_bs = candle_nn::Init::Uniform {
        lo: -bound,
        up: bound,
    };
    let bs = vb.get_with_hints(out_channels, "bias", init_bs)?;
    Ok(Conv2d::new(ws, Some(bs), cfg))
}

pub struct CNN {
    conv1: Conv2d,
    // act: Activation,
    // dropout: Dropout,
    conv2: Conv2d,
    profile: bool,
}

impl CNN {
    pub fn load(vb: VarBuilder, profile: bool) -> Result<Self> {
        let config = Conv2dConfig {
            padding: 0,
            stride: 1,
            dilation: 1,
            groups: 1,
        };

        let conv1 = conv2d_non_square(1024, 32, 7, 1, config, vb.pp("classifier.0"))?;
        // let act = Activation::Relu;
        // let dropout = Dropout::new(0.0);
        let conv2 = conv2d_non_square(32, 20, 7, 1, config, vb.pp("classifier.3"))?;
        // Ok(Self { conv1, act, dropout, conv2 })
        Ok(Self {
            conv1,
            conv2,
            profile,
        })
    }

    pub fn forward(&self, xs: &Tensor) -> Result<Tensor> {
        //println!("input shape: {:?}", xs.shape());
        let xs: Tensor = xs.permute((0, 2, 1))?.unsqueeze(D::Minus1)?;
        //println!("input after permutation: ");
        //println!("{xs}");
        //println!("{:?}", xs.shape());
        // println!("{:?}", xs.shape());
        let xs: Tensor = xs.pad_with_zeros(2, 3, 3)?;
        // println!("{:?}", xs.shape());
        let xs: Tensor = self.conv1.forward(&xs)?;
        // println!("{:?}", xs.shape());
        // println!("{xs}");
        //println!("{:?}", xs.shape());
        let xs: Tensor = xs.relu()?;
        let xs: Tensor = xs.clone();
        let xs: Tensor = xs.pad_with_zeros(2, 3, 3)?;
        let xs = self.conv2.forward(&xs)?.squeeze(D::Minus1)?;
        if self.profile {
            let xs: Tensor = xs.permute((0, 2, 1))?;
            let xs = log_softmax(&xs, D::Minus1)?;
            Ok(xs)
        } else {
            Ok(xs)
        }
    }
}
