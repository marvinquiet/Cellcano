import unittest

class TestImport(unittest.TestCase):
    def test_import_version(self):
        import Cellcano
        self.assertTrue(isinstance(Cellcano.__version__, str))
        self.assertNotEqual(len(Cellcano.__version__), 0)

    def test_import_submodule(self):
        from Cellcano.utils import _utils
        self.assertTrue(callable(_utils._csv_data_loader))

    def test_import_submodule_model(self):
        from Cellcano.models.distiller import Distiller
        self.assertTrue(callable(Distiller(None, None)))


if __name__ == '__main__':
    unittest.main()
