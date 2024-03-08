from inputCFALibrary import InputCFALibrary

if __name__ == '__main__':
    encoder = InputCFALibrary(ancillar_enhanced=False)
    # encoder.encode_on_the_fly()
    encoder.verify_horizontal_compatibility_for_reused_interface_cfas()