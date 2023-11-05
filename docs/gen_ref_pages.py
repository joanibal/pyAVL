"""Generate the code reference pages."""

from pathlib import Path

import mkdocs_gen_files

nav = mkdocs_gen_files.Nav()

src_dir = "moppy"
for path in sorted(Path(src_dir).rglob("*.py")):  #
    module_path = path.relative_to(src_dir).with_suffix("")  #
    if module_path.name == "__init__":
        continue

    doc_path = path.relative_to(src_dir).with_suffix(".md")  #
    full_doc_path = Path("reference", doc_path)  #

    parts = list(module_path.parts)

    if parts[-1] == "__init__":  #
        parts = parts[:-1]
    elif parts[-1] == "__main__":
        continue

    nav[parts] = doc_path.as_posix()  #

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:  #
        identifier = ".".join([src_dir] + parts)  #
        print("::: " + identifier, file=fd)  #

    mkdocs_gen_files.set_edit_path(full_doc_path, path)  #


with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:  #
    nav_file.writelines(nav.build_literate_nav())  #
