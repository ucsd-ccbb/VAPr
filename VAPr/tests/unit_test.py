import unittest

class Test_Utilities(unittest.TestCase):


    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

    def test_greater(self):
        argFoo = 123
        argBar = 452

        self.assertGreater(argFoo, argBar, "Foo is less than Bar")
        # -- this assert will fail

        self.assertLess(argFoo, argBar, "Foo is greater than Bar")
        # -- this assert will succeed

# Run the test case
if __name__ == '__main__':
    fooSuite = unittest.TestLoader().loadTestsFromTestCase(Test_Utilities)

    fooRunner = unittest.TextTestRunner()
    fooResult = fooRunner.run(fooSuite)

    print("---- START OF TEST RESULTS")
    print(fooResult)

    print("fooResult::errors")
    print(fooResult.errors)

    print("fooResult::failures")
    print(fooResult.failures)

    print("fooResult::skipped")
    print(fooResult.skipped)

    print("fooResult::successful")
    print(fooResult.wasSuccessful())

    print("fooResult::test-run")
    print(fooResult.testsRun)
    print("---- END OF TEST RESULTS")

