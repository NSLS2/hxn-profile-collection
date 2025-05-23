import copy

from event_model import (
    DocumentNames,
    schema_validators,
    unpack_datum_page,
    unpack_event_page,
)
from event_model.documents import (
    Datum,
    DatumPage,
    Event,
    EventDescriptor,
    EventPage,
    Resource,
    RunStart,
    RunStop,
    StreamDatum,
    StreamResource,
)

from bluesky.run_engine import Dispatcher
from bluesky.callbacks.core import CallbackBase
from bluesky.callbacks.tiled_writer import TiledWriter


def patch_descriptor(doc):
    # Issue #1: data keys for 'scaler_alive', 'sclr1_ch2' .. 'sclr1_ch8' have incorrect "dtype": "array"  and "shape": [1].
    #   Replace the values with "dtype": "number" and "shape": []. Also set "dtype_str": "<i8" for "sclr1_ch.." to get rid of warnings.
    # Issue #2: Fix for incorrect shape and type of merlin1 data.
    DATA_KEY_PREFIXES = ["scaler_alive", "sclr1_ch"]
    for k in doc["data_keys"].keys():
        if k == "scaler_alive":
            doc["data_keys"][k]["dtype"] = "number"
            doc["data_keys"][k]["shape"] = []
        elif k.startswith("sclr"):
            doc["data_keys"][k]["dtype"] = "number"
            doc["data_keys"][k]["dtype_str"] = "<f4"
        elif k == "merlin1":
            doc["data_keys"][k]["dtype_str"] = "<u4"
            # doc["data_keys"][k]["shape"] = [1, 209, 203]     # The shape is inconsistent; need to validate
        elif k == "eiger1" or k == "eiger_mobile_image":
            # doc["data_keys"][k]["shape"] = [1, 409, 412]     # or [1, 394, 338]? in "TPX_HDF5" spec
            doc["data_keys"][k]["dtype_str"] = "<u2"           # Inconsistent? Need to validate
        elif k == "eiger2_image":
            doc["data_keys"][k]["dtype_str"] = "<u4"           # Could be "<u4" or "<u2"
        elif k.startswith("Det"):
            doc["data_keys"][k]["dtype_str"] = "<f4"

    return doc

def patch_datum(doc):
    if "datum_kwargs" in doc.keys():
        dataset = doc['datum_kwargs'].pop("det_elem", None) \
                or doc['datum_kwargs'].pop("column", None) \
                or doc['datum_kwargs'].pop("field", None)
        if dataset:
            doc["datum_kwargs"]["dataset"] = dataset
        
        if channel := doc["datum_kwargs"].get("channel", None):
            doc["datum_kwargs"]["slice"] = f"(:,{channel-1},:)"
            doc["datum_kwargs"]["squeeze"] = True
    return doc

def patch_resource(doc):
    if doc["resource_path"].startswith(doc["root"]):
        doc["resource_path"] = doc["resource_path"].replace(doc["root"], '')
    # doc["root"] == re.sub(r"^/data", "/nsls2/data/hxn/legacy", doc["root"])
    doc["root"] = doc["root"].replace("/data", "/nsls2/data/hxn/legacy")
    doc["resource_path"] = doc["resource_path"].replace("nsls2/data2/hxn", "nsls2/data/hxn")
    
    # Replace 'frame_per_point' with 'multiplier'
    if frame_per_point := doc["resource_kwargs"].pop("frame_per_point", None):
        doc["resource_kwargs"]["multiplier"] = frame_per_point
    
    # Sey default chunk_shape to (1,)
    if "chunk_shape" not in doc["resource_kwargs"].keys():
        doc["resource_kwargs"]["chunk_shape"] = (1, )

    spec = doc["spec"]
    if spec in ["MERLIN_FLY_STREAM_V2", "TPX3", "TPX_HDF5", "EIGER2_STREAM", "MERLIN_HDF5_BULK"]:
        doc["resource_kwargs"].update({"dataset": '/entry/data/data', "chunk_shape": (1, ), "join_method": "concat"})
    elif spec == "PANDA":
        doc["resource_kwargs"].update({"chunk_shape": (1024, ), "join_method": "concat"})
    elif spec == "ROI_HDF5_FLY":
        doc["resource_kwargs"].update({"chunk_shape": (1024, ), "join_method": "concat"})
    elif spec in ["XSP"]:
        # NOTE: here data is accessed by "channel"
        doc["resource_kwargs"].update({"dataset": 'entry/instrument/detector/data', "chunk_shape": (1, ), "join_method": "concat"})
    elif spec == "XSP3_BULK":
        # NOTE: here data is accessed by "channel"
        # NOTE: "XSP3_BULK" does not declare the number of rows in the descriptor -- need to validate (see below)
        doc["resource_kwargs"].update({"dataset": 'entry/instrument/detector/data',
                                       "chunk_shape": (1, ),
                                       "join_method": "stack"})
    elif spec == "XSP3":
        # NOTE: here data is accessed by "channel" -- need to decide how to handle it in Tiled
        doc["resource_kwargs"].update({"dataset": 'entry/instrument/detector/data', "chunk_shape": (1, ), "join_method": "stack"})
    elif spec == "SIS_HDF51_FLY_STREAM_V1":
        # These correspond to 'sclr1_ch*' data_keys
        doc["resource_kwargs"].update({"chunk_shape": (1024, ), "join_method": "concat"})

    # Validate the structure if e.g. the datum shape is not declared in the descriptor
    # NOTE: Flyscannning specs often do not include the total number recorded points, only points per frame
    # NOTE: TPX_HDF5 has inconsistent shape definitions for eiger1 and merlin (stackable T/F)
    if spec in ["XSP3_BULK", "MERLIN_HDF5_BULK", "ROI_HDF5_FLY", "SIS_HDF51_FLY_STREAM_V1",
                "MERLIN_FLY_STREAM_V2", "PANDA"] + ["TPX_HDF5"]:
        doc["resource_kwargs"].update({"_validate": True})

    return doc


class DocumentConverter(CallbackBase):
    """Callback for updating Bluesky documents"""

    def __init__(self):
        self._token_refs = {}
        self.dispatcher = Dispatcher()

    def start(self, doc: RunStart):
        doc = copy.deepcopy(doc)
        self.emit(DocumentNames.start, doc)

    def stop(self, doc: RunStop):
        doc = copy.deepcopy(doc)
        self.emit(DocumentNames.stop, doc)

    def descriptor(self, doc: EventDescriptor):
        doc = copy.deepcopy(doc)
        doc = patch_descriptor(doc)
        self.emit(DocumentNames.descriptor, doc)

    def event(self, doc: Event):
        doc = copy.deepcopy(doc)
        self.emit(DocumentNames.event, doc)

    def resource(self, doc: Resource):
        doc = copy.deepcopy(doc)
        doc = patch_resource(doc)
        self.emit(DocumentNames.resource, doc)

    def stream_resource(self, doc: StreamResource):
        doc = copy.deepcopy(doc)
        self.emit(DocumentNames.stream_resource, doc)

    def datum(self, doc: Datum):
        doc = copy.deepcopy(doc)
        doc = patch_datum(doc)
        self.emit(DocumentNames.datum, doc)

    def stream_datum(self, doc: StreamDatum):
        doc = copy.deepcopy(doc)
        self.emit(DocumentNames.stream_datum, doc)

    def datum_page(self, doc: DatumPage):
        for _doc in unpack_datum_page(doc):
            self.datum(_doc)

    def event_page(self, doc: EventPage):
        for _doc in unpack_event_page(doc):
            self.event(_doc)

    def emit(self, name, doc):
        """Check the document schema and send to the dispatcher"""
        schema_validators[name].validate(doc)
        self.dispatcher.process(name, doc)

    def subscribe(self, func, name="all"):
        """Convenience function for dispatcher subscription"""
        token = self.dispatcher.subscribe(func, name)
        self._token_refs[token] = func
        return token

    def unsubscribe(self, token):
        """Convenience function for dispatcher un-subscription"""
        self._token_refs.pop(token, None)
        self.dispatcher.unsubscribe(token)


converter = DocumentConverter()
tw = TiledWriter(client= ...)
converter.subscribe(tw)

RE.subscribe(converter)