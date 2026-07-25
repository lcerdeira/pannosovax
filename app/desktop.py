#!/usr/bin/env python3
"""PanNosoVax Studio — lançador desktop (janela nativa via pywebview).

Sobe o backend FastAPI numa thread e abre uma janela nativa apontando para ele.
É o ponto de entrada empacotado (`pannosovax-studio`); no macOS/Windows o usuário
abre um app, sem terminal.

Dev:  python app/desktop.py
"""
from __future__ import annotations
import socket
import threading
import time

HOST = "127.0.0.1"


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((HOST, 0))
        return s.getsockname()[1]


def _serve(port: int):
    import uvicorn
    from app.backend import app  # import tardio: só quando a janela vai abrir
    uvicorn.run(app, host=HOST, port=port, log_level="warning")


def _wait_up(port: int, timeout: float = 20.0) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection((HOST, port), timeout=0.5):
                return True
        except OSError:
            time.sleep(0.2)
    return False


def main():
    import webview  # pywebview

    port = _free_port()
    threading.Thread(target=_serve, args=(port,), daemon=True).start()
    if not _wait_up(port):
        raise SystemExit("backend não subiu a tempo")
    webview.create_window("PanNosoVax Studio", f"http://{HOST}:{port}",
                          width=1024, height=760, min_size=(800, 600))
    webview.start()


if __name__ == "__main__":
    main()
