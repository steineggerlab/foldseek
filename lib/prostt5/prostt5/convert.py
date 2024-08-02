import ctranslate2
try:
    converter = ctranslate2.converters.TransformersConverter("Rostlab/ProstT5_fp16", load_as_float16=True)
except Exception as e:
    print(f"Failed to load converter: {e}")

converter.convert("ProstT5_ct2_fp16", force=True)
