import io
import tarfile
import tempfile
from pathlib import Path
import unittest
from prepare_taxonomic_lineages import ancestry, lowest_common_ancestor, read_taxdump, resolve_taxid


class TaxonomyTests(unittest.TestCase):
    def test_parsing_resolution_and_lca(self):
        files = {
            "nodes.dmp": "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tdomain\t|\n3\t|\t2\t|\tfamily\t|\n4\t|\t3\t|\tgenus\t|\n5\t|\t3\t|\tgenus\t|\n",
            "names.dmp": "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tLife\t|\t\t|\tscientific name\t|\n",
            "merged.dmp": "9\t|\t4\t|\n", "delnodes.dmp": "8\t|\n"}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "tax.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                for name, text in files.items():
                    payload = text.encode(); info = tarfile.TarInfo(name); info.size = len(payload)
                    archive.addfile(info, io.BytesIO(payload))
            parents, ranks, names, merged, deleted = read_taxdump(path)
            self.assertEqual(resolve_taxid("9", parents, merged, deleted), "4")
            self.assertEqual(ancestry("4", parents), ["4", "3", "2", "1"])
            self.assertEqual(lowest_common_ancestor(["4", "5"], parents), "3")
            with self.assertRaises(ValueError): resolve_taxid("8", parents, merged, deleted)


if __name__ == "__main__": unittest.main()
