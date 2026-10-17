from typing import Dict

from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware

from app.api.endpoints import router as api_router
from app.api.ground_endpoints import router as ground_router

app = FastAPI(title="SSE_H2SQ Simulation API", version="1.0.0")

# Enable CORS for the frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Adjust in production to specific frontend domains
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(api_router, prefix="/api")
app.include_router(ground_router, prefix="/api")


@app.get("/")
def read_root():
    return {"message": "Welcome to the SSE_H2SQ Simulation API"}


# Global dictionary to keep track of active websocket connections
active_connections: Dict[str, WebSocket] = {}


@app.websocket("/ws/{client_id}")
async def websocket_endpoint(websocket: WebSocket, client_id: str):
    await websocket.accept()
    active_connections[client_id] = websocket
    try:
        while True:
            # Keep the connection alive and listen for
            # client messages if needed
            await websocket.receive_text()
    except WebSocketDisconnect:
        if client_id in active_connections:
            del active_connections[client_id]
