"""Check every equation label, displayed tag, reference link, and source body."""

from collections import Counter
from pathlib import Path
import difflib
import hashlib
import json
import re
import fitz
from pypdf import PdfReader

REPORT = Path(__file__).resolve().parent
ROOT = REPORT.parents[3]
SOURCE = ROOT / "research-candidates/vol2_boundary026/collision006"
OLD_SOURCE = ROOT / "research-candidates/vol2_boundary026/collision005"
PDF = REPORT / "build/reader.pdf"
OLD_PDF = ROOT / "reports/research/vol2_boundary026/collision005/build-pdftex/reader.pdf"


def groups(text):
    depth = 0
    result = []
    for i, char in enumerate(text):
        if char == "{":
            if depth == 0:
                start = i + 1
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                result.append(text[start:i])
    return result


def labels(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith(r"\newlabel{"):
            key, value = groups(line)
            if not key.endswith("@cref"):
                result[key] = groups(value)
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def link_rows(document, reader):
    rows = []
    for index, page in enumerate(reader.pages):
        for raw in page.get("/Annots", []):
            annotation = raw.get_object()
            action = annotation.get("/A", {})
            if action.get("/S") != "/GoTo":
                continue
            x0, y0, x1, y1 = map(float, annotation["/Rect"])
            height = document[index].rect.height
            area = fitz.Rect(x0 - .3, height - y1 - .3, x1 + .3, height - y0 + .3)
            visible = re.sub(r"\s+", "", document[index].get_text(clip=area))
            rows.append({"page": index + 1, "destination": str(action["/D"]),
                         "text": visible, "rect": list(area)})
    return rows


def visible_text_without_references(document, links, tags):
    excluded = {i: [] for i in range(len(document))}
    for link in links:
        # Only reference-number annotations are removed. Prose remains.
        if re.fullmatch(r"[\[\]().,\d]+", link["text"]):
            excluded[link["page"] - 1].append(fitz.Rect(link["rect"]))
    for tag in tags:
        excluded[tag["page"] - 1].append(fitz.Rect(tag["rect"]))
    output = []
    for index, page in enumerate(document):
        for block in page.get_text("rawdict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    for char in span["chars"]:
                        box = fitz.Rect(char["bbox"])
                        center = (box.tl + box.br) / 2
                        if center.y < 75 or center.y > 739:
                            continue
                        if any(rect.contains(center) for rect in excluded[index]):
                            continue
                        output.append(char["c"])
    return re.sub(r"\s+", "", "".join(output))


def run():
    body_checks = []
    source_refs = Counter()
    source_keys = []
    for path in sorted((SOURCE / "chapters").glob("*.tex")):
        predecessor = OLD_SOURCE / path.relative_to(SOURCE)
        assert path.read_bytes() == predecessor.read_bytes(), path
        body_checks.append({"path": str(path.relative_to(SOURCE)), "sha256": sha(path)})
        text = path.read_text()
        source_keys += re.findall(r"\\label\{(eq:[^}]+)\}", text)
        source_refs.update(re.findall(r"\\(?:eqref|refeq|ref)\{(eq:[^}]+)\}", text))
    assert len(source_keys) == len(set(source_keys)) == 30
    doc = fitz.open(PDF)
    olddoc = fitz.open(OLD_PDF)
    reader = PdfReader(PDF)
    oldreader = PdfReader(OLD_PDF)
    aux = labels(PDF.with_suffix(".aux"))
    links = link_rows(doc, reader)
    oldlinks = link_rows(olddoc, oldreader)
    assert len(doc) == len(olddoc) == 60
    equation_rows = []
    sequence = Counter()
    tag_regions = []
    for key in [key for key in aux if key in source_keys]:
        number, printed_page, title, destination, _ = aux[key]
        chapter, counter = map(int, number.split("."))
        sequence[chapter] += 1
        assert counter == sequence[chapter], (key, number, sequence)
        assert destination.startswith("equation."), (key, destination)
        target = reader.named_destinations[destination]
        index = reader.get_destination_page_number(target)
        assert int(printed_page) == index + 1
        assert str(target.typ) == "/XYZ"
        y = doc[index].rect.height - float(target.top)
        candidates = [word for word in doc[index].get_text("words")
                      if word[4] == f"({number})" and word[0] > .7 * doc[index].rect.width
                      and -2 < word[1] - y < 22]
        assert len(candidates) == 1, (key, candidates, y)
        tag = candidates[0]
        area = fitz.Rect(tag[:4])
        overlap = [word for word in doc[index].get_text("words")
                   if word != tag and fitz.Rect(word[:4]).intersects(area)]
        assert not overlap, (key, overlap)
        refs = [row for row in links if row["destination"] == destination]
        assert len(refs) == source_refs[key], (key, len(refs), source_refs[key])
        assert all(row["text"].strip("[]().,") == number for row in refs), (key, refs)
        tag_regions.append({"page": index + 1, "rect": list(area)})
        equation_rows.append({"label": key, "number": number,
                              "destination": destination, "page": index + 1,
                              "target_top": float(target.top), "tag_rect": list(area),
                              "reference_links": refs})
    assert len({row["destination"] for row in equation_rows}) == 30
    assert sequence == {1: 3, 2: 2, 3: 8, 4: 11, 5: 6}
    before = visible_text_without_references(olddoc, oldlinks, [])
    after = visible_text_without_references(doc, links, tag_regions)
    (REPORT / "text-before-normalized.txt").write_text(before)
    (REPORT / "text-after-normalized.txt").write_text(after)
    prefix = 0
    while prefix < min(len(before), len(after)) and before[prefix] == after[prefix]:
        prefix += 1
    suffix = 0
    while suffix < min(len(before), len(after)) - prefix and before[-1-suffix] == after[-1-suffix]:
        suffix += 1
    a_end, b_end = len(before) - suffix, len(after) - suffix
    matcher = difflib.SequenceMatcher(a=before[prefix:a_end], b=after[prefix:b_end], autojunk=False)
    changes = []
    for kind, i, j, k, l in matcher.get_opcodes():
        if kind == "equal":
            continue
        i, j, k, l = i + prefix, j + prefix, k + prefix, l + prefix
        changes.append({"kind": kind, "before": before[i:j], "after": after[k:l],
                        "before_context": before[max(0, i-70):min(len(before), j+70)],
                        "after_context": after[max(0, k-70):min(len(after), l+70)]})
    assert all(row["kind"] == "insert" and row["after"] == "-" for row in changes)
    line_break_checks = []
    for page_number, ending, continuation in [(53, "correspond-", "ing"), (55, "pre-", "cisely"), (56, "con-", "traction")]:
        lines = []
        for block in doc[page_number - 1].get_text("rawdict")["blocks"]:
            for line in block.get("lines", []):
                lines.append("".join(char["c"] for span in line["spans"] for char in span["chars"]))
        hits = [n for n, line in enumerate(lines[:-1]) if line.endswith(ending) and lines[n+1].startswith(continuation)]
        assert len(hits) == 1, (page_number, ending, hits)
        line_break_checks.append({"page": page_number, "line_end": lines[hits[0]], "next_line": lines[hits[0]+1]})
    assert len(changes) == len(line_break_checks) == 3
    result = {"pdf_sha256": sha(PDF), "pages": len(doc), "unchanged_mathematical_bodies": body_checks,
              "equations": equation_rows, "equation_count": len(equation_rows),
              "equation_reference_link_count": sum(len(row["reference_links"]) for row in equation_rows),
              "number_sequence_by_chapter": dict(sequence),
              "proof_text_equal_after_removing_numeric_references_tags_headers": before == after,
              "remaining_text_differences": changes,
              "line_break_hyphenation_checks": line_break_checks,
              "proof_text_preserved_modulo_verified_reference_fields_and_three_line_break_hyphens": True,
              "extraction_limit": "Whitespace, numeric internal-reference rectangles, equation tags, and page header/footer bands are excluded. Source byte identity is separately checked for all eight mathematical bodies."}
    (REPORT / "reader-checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ["pages", "equation_count", "equation_reference_link_count", "number_sequence_by_chapter", "proof_text_equal_after_removing_numeric_references_tags_headers", "remaining_text_differences"]}, indent=2))


if __name__ == "__main__":
    run()
