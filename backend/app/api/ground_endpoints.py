# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
import json
import asyncio
from fastapi import APIRouter,HTTPException
from fastapi.responses import FileResponse,Response
from starlette.concurrency import run_in_threadpool
from app.services import ground_validation as service

router=APIRouter(prefix='/ground-validation')
_figure_lock=asyncio.Lock()


@router.get('')
async def studies():return await run_in_threadpool(service.list_jobs)


@router.post('')
async def create(request:service.Request):
    try:return await service.start(request)
    except (RuntimeError,OSError) as error:raise HTTPException(409,str(error)) from error


@router.get('/{identity}')
async def study(identity:str):
    try:return await run_in_threadpool(service.get,identity)
    except (FileNotFoundError,ValueError) as error:raise HTTPException(404,'Study not found') from error


@router.post('/{identity}/cancel')
async def cancel(identity:str):return dict(cancelled=service.cancel(identity))


@router.post('/{identity}/resume')
async def resume(identity:str):
    try:return await service.resume(identity)
    except (RuntimeError,OSError) as error:raise HTTPException(409,str(error)) from error


@router.get('/{identity}/export')
async def export(identity:str,format:str='csv',point:str='',observable:str='E',quantity:str='mean'):
    job=await study(identity)
    if format=='json':return Response(json.dumps(job,indent=2,allow_nan=False),media_type='application/json')
    if format in ('pdf','svg','png') and point in {p['id'] for p in job['points'] if p['levels']} and observable in service.ground_ed.KEYS and quantity in ('mean','difference','sampling_difference','precision','tau_int','ess','rhat','fidelity','state_infidelity','state_tvd','state_missing'):
        async with _figure_lock:
            try:content=await run_in_threadpool(service.figure_result,job,point,observable,quantity,format)
            except ValueError as error:raise HTTPException(422,str(error)) from error
        return Response(content,media_type={'pdf':'application/pdf','svg':'image/svg+xml','png':'image/png'}[format])
    if format!='csv':raise HTTPException(422,'Choose csv or json')
    return Response(service.csv_result(job),media_type='text/csv',headers={'Content-Disposition':f'attachment; filename="ground-validation-{identity[:8]}.csv"'})


@router.get('/{identity}/traces/{point}/{projection}/{method}/{chain}')
async def trace(identity:str,point:str,projection:float,method:str,chain:int):
    job=await study(identity)
    try:path=service.trace_path(job,point,projection,method,chain)
    except (FileNotFoundError,ValueError) as error:raise HTTPException(404,'Trace not found') from error
    if not path.is_file():raise HTTPException(404,'Trace not ready')
    return FileResponse(path,filename=f'ground-{point}-{projection:g}-{method}-{chain}.npz')
