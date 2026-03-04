"""
Standalone timmNet module for Open-Canopy inference.

Embeds the minimal classes from the Open-Canopy repository so that PVTv2
inference works without cloning the external repo.

Original source: https://github.com/fajwel/Open-Canopy
  - src/models/components/timmNet.py
  - src/models/components/utils/utils.py   (infer_output, set_first_layer)
  - src/models/components/utils/seg_blocks.py (SimpleSegmentationHead)
"""

import math

import timm
import torch
import torch.nn as nn
import torch.nn.functional as F


# ======================================================================
# From src/models/components/utils/utils.py
# ======================================================================

def infer_output(model, num_channels, img_size, **kwargs):
    """Infer the output dimensions of the backbone model."""
    dummy_image = torch.rand(
        1, num_channels, img_size, img_size, dtype=torch.float32
    )
    with torch.no_grad():
        dummy_features = model.forward_features(dummy_image, **kwargs)

    def analyse(img_size, dummy_features):
        remove_cls_token = False
        print(f"backbone output features.shape: {dummy_features.shape}")
        if len(dummy_features.shape) == 4:
            embed_dim = dummy_features.shape[-3]
            downsample_factor = img_size // dummy_features.shape[-1]
            feature_size = dummy_features.shape[-1]
            features_format = "NCHW"
        elif len(dummy_features.shape) == 3:
            embed_dim = dummy_features.shape[-1]
            feature_size = math.floor(math.sqrt(dummy_features.shape[-2]))
            if feature_size ** 2 != dummy_features.shape[-2]:
                if feature_size ** 2 == dummy_features.shape[-2] - 1:
                    remove_cls_token = True
                else:
                    raise ValueError(
                        f"backbone output features.shape[-2] must be a square "
                        f"number? Currently it is {dummy_features.shape[-2]}"
                    )
            downsample_factor = img_size // feature_size
            features_format = "NLC"
        else:
            raise ValueError(
                f"backbone output features.shape must be of dimension 3 or 4? "
                f"Currently it is {dummy_features.shape}"
            )
        return (embed_dim, downsample_factor, feature_size,
                features_format, remove_cls_token)

    if isinstance(dummy_features, (list, tuple)):
        (embed_dim, downsample_factor, feature_size,
         features_format, remove_cls_token) = zip(
            *[analyse(img_size, f) for f in dummy_features]
        )
        assert all(cls == remove_cls_token[0] for cls in remove_cls_token)
        assert all(f == features_format[0] for f in features_format)
        embed_dim = list(embed_dim)
        downsample_factor = list(downsample_factor)
        feature_size = list(feature_size)
        features_format = features_format[0]
        remove_cls_token = remove_cls_token[0]
    elif isinstance(dummy_features, torch.Tensor):
        (embed_dim, downsample_factor, feature_size,
         features_format, remove_cls_token) = analyse(img_size, dummy_features)

    return (embed_dim, downsample_factor, feature_size,
            features_format, remove_cls_token)


def set_first_layer(model, n_channels, is_rgb=None):
    """Adapt the first layer of a model from 3 channels to n_channels."""
    if n_channels == 3:
        return

    if is_rgb is None:
        is_rgb = n_channels >= 3

    if is_rgb:
        assert n_channels > 3

    for module in model.modules():
        if isinstance(module, nn.Conv2d) and module.in_channels == 3:
            break
        if isinstance(module, nn.Linear):
            break
    previous_weight = module.weight.detach()

    if previous_weight.dim() == 4:
        n_out = previous_weight.shape[0]
        if is_rgb:
            new_weight = torch.randn(
                n_out, n_channels,
                previous_weight.shape[2], previous_weight.shape[3]
            )
            new_weight[:, :3] = previous_weight
        else:
            mean = previous_weight.mean(dim=1)
            new_weight = torch.stack([mean] * n_channels, dim=1)
    elif previous_weight.dim() == 2:
        n_out = previous_weight.shape[0]
        n_elem = previous_weight.shape[1] // 3
        if is_rgb:
            new_weight = torch.randn((n_out, n_channels * n_elem))
            new_weight[:, :3 * n_elem] = previous_weight
        else:
            mean = previous_weight.reshape(n_out, -1, 3).mean(dim=-1)
            new_weight = torch.stack([mean] * n_channels, dim=-1)
            new_weight = new_weight.reshape(n_out, -1)

    module.weight = nn.parameter.Parameter(new_weight)


# ======================================================================
# From src/models/components/utils/seg_blocks.py
# ======================================================================

class SimpleSegmentationHead(nn.Module):
    """Simple segmentation head with transposed convolutions."""

    def __init__(
        self,
        embed_dim,
        downsample_factor,
        remove_cls_token,
        features_format,
        features_sizes,
        num_classes,
        decoder_stride=2,
        **kwargs,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.downsample_factor = downsample_factor
        self.remove_cls_token = remove_cls_token
        self.features_format = features_format
        self.feature_size = features_sizes
        self.num_classes = num_classes
        self.decoder_stride = decoder_stride

        self.layered_output = isinstance(self.embed_dim, (list, tuple))
        if self.layered_output:
            self.embed_dim = self.embed_dim[-1]
            self.downsample_factor = self.downsample_factor[-1]
            self.feature_size = self.feature_size[-1]

        print(f"{self.embed_dim=}, {self.downsample_factor=}, "
              f"{self.feature_size=}")

        depth = math.log(self.downsample_factor, decoder_stride)
        assert depth.is_integer(), (
            f"decoder stride({decoder_stride}) must be a power of the "
            f"downsample factor({self.downsample_factor})"
        )
        depth = int(depth)

        self.layers = nn.Sequential(
            *[
                nn.Sequential(
                    nn.ConvTranspose2d(
                        self.embed_dim // 2 ** d,
                        self.embed_dim // 2 ** (d + 1),
                        decoder_stride,
                        stride=decoder_stride,
                    ),
                    nn.BatchNorm2d(self.embed_dim // 2 ** (d + 1)),
                    nn.GELU(),
                    nn.Conv2d(
                        self.embed_dim // 2 ** (d + 1),
                        self.embed_dim // 2 ** (d + 1),
                        3,
                        padding="same",
                    ),
                    nn.BatchNorm2d(self.embed_dim // 2 ** (d + 1)),
                    nn.GELU(),
                )
                for d in range(depth - 1)
            ]
            + [
                nn.ConvTranspose2d(
                    self.embed_dim // 2 ** (depth - 1),
                    num_classes,
                    decoder_stride,
                    stride=decoder_stride,
                )
            ]
        )

    def forward(self, x):
        if self.layered_output:
            x = x[-1]
        if self.remove_cls_token:
            x = x[:, 1:, :]
        if self.features_format == "NLC":
            x = x.reshape(
                x.shape[0], self.feature_size, self.feature_size, x.shape[-1]
            )
            x = x.permute(0, 3, 1, 2)
        return self.layers(x)


# ======================================================================
# From src/models/components/timmNet.py
# ======================================================================

class timmNet(nn.Module):
    """Segmentation network using timm models as backbone.

    Standalone version that does not require the Open-Canopy source tree.
    """

    def __init__(
        self,
        backbone="vit_base_patch16_384",
        num_classes=4,
        num_channels=1,
        segmentation_head=SimpleSegmentationHead,
        pretrained=True,
        pretrained_path=None,
        img_size=512,
        lora_rank=0,
        use_FPN=False,
        chkpt_path=None,
        decoder_stride=None,
    ):
        super().__init__()
        self.backbone = backbone
        self.num_classes = num_classes
        self.num_channels = num_channels
        self.pretrained = pretrained
        self.pretrained_path = pretrained_path
        self.img_size = img_size

        if pretrained_path is not None:
            pretrained_cfg_overlay = dict(file=pretrained_path)
        else:
            pretrained_cfg_overlay = dict()

        print(f"Using timm version: {timm.__version__}")

        if backbone.startswith("swin"):
            additional_arg = {
                "img_size": self.img_size,
                "features_only": True,
            }
        elif backbone.startswith("pvt_v2"):
            additional_arg = {"features_only": True}
        elif backbone.startswith("twins_pcpvt"):
            additional_arg = {}
        elif backbone.startswith("vit_base_r50"):
            additional_arg = {"num_classes": 0}
        else:
            additional_arg = {}

        self.model = timm.create_model(
            self.backbone,
            pretrained=self.pretrained,
            in_chans=3,
            pretrained_cfg_overlay=pretrained_cfg_overlay,
            **additional_arg,
        )

        set_first_layer(self.model, num_channels)

        if not hasattr(self.model, "forward_features"):
            self.model.forward_features = self.model.forward

        if backbone.startswith("swin"):
            class reorder_swin(nn.Module):
                def __init__(self, model):
                    super().__init__()
                    self.model = model

                def forward(self, x):
                    return [
                        feature.permute(0, 3, 1, 2)
                        for feature in self.model(x)
                    ]

                def forward_features(self, x):
                    return self.forward(x)

            self.model = reorder_swin(self.model)

        (
            self.embed_dim,
            self.downsample_factor,
            self.feature_size,
            self.features_format,
            self.remove_cls_token,
        ) = infer_output(self.model, self.num_channels, self.img_size)

        # decoder_stride : si None, utiliser downsample_factor (1 seule couche)
        _ds = decoder_stride if decoder_stride is not None else self.downsample_factor
        self.seg_head = segmentation_head(
            self.embed_dim,
            self.downsample_factor,
            self.remove_cls_token,
            self.features_format,
            self.feature_size,
            self.num_classes,
            decoder_stride=_ds,
        )

    def forward(self, x, metas=None):
        x = self.model.forward_features(x)
        x = self.seg_head(x)
        return {"out": x}
