"""Écriture TOML : aller-retour exact sur tous les exemples."""


def test_toml_roundtrip_examples():
    from microrans.cli import examples_dir
    from microrans.mesh2d.builder import load_config
    from microrans.tomlio import dumps, loads
    for f in examples_dir().glob("*.toml"):
        cfg = load_config(f)
        assert loads(dumps(cfg)) == cfg, f.name
