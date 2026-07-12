# Import from the submodules directly; keeping this __init__ empty means importing a constant
# never drags in the heavy analysis stack. This lets the Pixi eager env, which must use an
# old version of nextflow, import evaluation_utils.constants without error (see run_eager.py).
