import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "chapters/theory/sc_chtop_heptagon.tex"


def read() -> str:
    return TARGET.read_text()


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def test_sc_chtop_definition_is_chain_coloured_dioperad_not_only_name():
    body = compact(read())

    required = (
        r"C_\bullet^{\log}\SCchtop(\mathsf{cl}^{p},\mathsf{op}^{q};\mathsf{cl})",
        r"C_\bullet^{\log}\SCchtop(\mathsf{cl}^{p},\mathsf{op}^{q};\mathsf{op})",
        r"C_\bullet^{\log}\bigl(\FM_p(\C)\times\Conf_q^{<}(\R)\bigr)",
        r"theopen-to-closedcomponentisthezerochaincomplex",
    )
    for needle in required:
        assert needle in body


def test_sc_chtop_definition_prints_local_composition_maps():
    body = compact(read())

    required = (
        r"\gamma_{\mathsf{cl}}\colon\FM_p(\C)\times\prod_{i=1}^{p}\FM_{n_i}(\C)",
        r"\gamma_{\mathsf{op}}\colon\bigl(\FM_p(\C)\times\Conf_q^{<}(\R)\bigr)",
        r"\times\prod_{j=1}^{q}\Conf_{m_j}^{<}(\R)",
        r"\FM_{\sum_in_i}(\C)\times\Conf_{\sum_jm_j}^{<}(\R)",
    )
    for needle in required:
        assert needle in body


def test_sc_chtop_w_resolution_has_full_orientation_line_and_differential():
    body = compact(read())

    required = (
        r"\mathfrak{o}(T)=\det\bigl(\mathbb{k}^{E_{\mathrm{int}}(T)}\bigr)\otimes\bigotimes_{v\inV(T)}\mathfrak{o}(v)",
        r"C_\bullet\!\left(W(\SCchtop)(\mathbfa;b)\right)",
        r"\bigoplus_{T:\mathbfa\tob}",
        r"C_\bullet\bigl([0,1]^{E_{\mathrm{int}}(T)}\bigr)",
        r"d_W=d_{\mathrm{chain}}+d_{[0,1]}",
        r"d_{\mathrm{contract},e}",
        r"d_{\mathrm{compose},e}",
        r"vertex-orientationtensorfactors",
    )
    for needle in required:
        assert needle in body


def test_sc_chtop_definition_rejects_unoriented_w_resolution():
    text = read()

    retired = (
        "tree enumeration without orientation",
        "orientation line optional",
        "some coherence",
    )
    for phrase in retired:
        assert phrase not in text
