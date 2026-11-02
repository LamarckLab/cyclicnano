"""Read-only web browser for cyclicnano design runs.

Importing this package pulls in nothing: `store` is standard library only and
`app` needs FastAPI, so each is imported on its own by whoever needs it. That
way the reading layer stays testable in an environment without a web stack.

Run it with `python -m cyclicnano.webui --root /data/lmk/cyclicnano_outputs`.
"""
