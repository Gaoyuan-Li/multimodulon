#!/usr/bin/env python3
"""
Sanity check for multimodulon imports and PyTorch dependencies.
"""

def test_pytorch_imports():
    """Test if PyTorch and related dependencies can be imported."""
    try:
        import torch
        print(f"✓ PyTorch {torch.__version__} imported successfully")
        print(f"  CUDA available: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"  CUDA device: {torch.cuda.get_device_name(0)}")
    except ImportError as e:
        print(f"✗ PyTorch import failed: {e}")
        return False
    
    try:
        import geotorch
        print(f"✓ geotorch imported successfully")
    except ImportError as e:
        print(f"✗ geotorch import failed: {e}")
        return False
    
    try:
        import scipy.optimize
        print(f"✓ scipy.optimize imported successfully")
    except ImportError as e:
        print(f"✗ scipy.optimize import failed: {e}")
        return False
    
    return True


def test_core_dependencies():
    """Test if core dependencies can be imported."""
    deps = [
        ('pandas', 'pd'),
        ('numpy', 'np'),
        ('matplotlib.pyplot', 'plt'),
        ('sklearn.model_selection', None),
    ]
    
    all_good = True
    for dep, alias in deps:
        try:
            if alias:
                exec(f"import {dep} as {alias}")
            else:
                exec(f"import {dep}")
            print(f"✓ {dep} imported successfully")
        except ImportError as e:
            print(f"✗ {dep} import failed: {e}")
            all_good = False
    
    return all_good


def test_multimodulon_imports():
    """Test if multimodulon modules can be imported."""
    try:
        from multimodulon.multiview_ica import run_multiview_ica_native
        print("✓ run_multiview_ica_native imported successfully")
    except ImportError as e:
        print(f"✗ run_multiview_ica_native import failed: {e}")
        return False
    
    try:
        from multimodulon.multiview_ica_optimization import run_nre_optimization_native
        print("✓ run_nre_optimization_native imported successfully")
    except ImportError as e:
        print(f"✗ run_nre_optimization_native import failed: {e}")
        return False
    
    try:
        from multimodulon import MultiModulon
        print("✓ MultiModulon class imported successfully")
    except ImportError as e:
        print(f"✗ MultiModulon class import failed: {e}")
        return False
    
    return True


def main():
    print("="*60)
    print("MULTIMODULON DEPENDENCY SANITY CHECK")
    print("="*60)
    
    print("\n1. Testing PyTorch and related dependencies...")
    pytorch_ok = test_pytorch_imports()
    
    print("\n2. Testing core scientific dependencies...")
    core_ok = test_core_dependencies()
    
    print("\n3. Testing multimodulon imports...")
    multimodulon_ok = test_multimodulon_imports()
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    if pytorch_ok and core_ok and multimodulon_ok:
        print("✓ All dependencies are available!")
        print("✓ The multimodulon package should work correctly.")
        return True
    else:
        print("✗ Some dependencies are missing.")
        if not pytorch_ok:
            print("  Install PyTorch: pip install torch==2.6.0 torchvision==0.21.0 torchaudio==2.6.0 geotorch==0.3.0 --index-url https://download.pytorch.org/whl/cu124")
        if not core_ok:
            print("  Install core dependencies: pip install pandas numpy matplotlib scikit-learn scipy")
        return False


if __name__ == "__main__":
    main()