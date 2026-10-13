"""Compatibility application factory for the FastAPI service.

The authoritative ASGI application is agents.api:app.
"""


def create_app():
    try:
        from agents.api import app
    except ModuleNotFoundError as exc:
        if exc.name in {"fastapi", "pydantic"}:
            return None
        raise
    return app
