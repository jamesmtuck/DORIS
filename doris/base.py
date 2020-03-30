from doris.base_conversion import convertIntToBytes,convertBytesToInt

class BaseCodec(object):
    def __init__(self,CodecObj=None,Policy=None):
        self._Obj = CodecObj
        if Policy is None:
            self._Policy = AllowAll()
        else:
            self._Policy = Policy
    def _encode(self, s):
        return s
    def encode(self,s):
        if self._Obj != None:
            return self._encode(self._Obj.encode(s))
        else:
            return self._encode(s)
    def _decode(self, s):
        return s
    def decode(self,s):
        s = self._decode(s)
        if self._Obj != None:
            return self._Obj.decode(s)
        else:
            return s

class TableCodec(BaseCodec):
    def __init__(self,CodecObj=None,keyEncWidth=20,keyDecWidth=4,cwEncWidth=5,cwDecWidth=1,Policy=None):
        BaseCodec.__init__(self,CodecObj=CodecObj,Policy=Policy)
        self._keyEncWidth = keyEncWidth
        self._keyDecWidth = keyDecWidth
        self._cwEncWidth = cwEncWidth
        self._cwDecWidth = cwDecWidth
        assert keyEncWidth % cwEncWidth == 0
        assert keyEncWidth / cwEncWidth == keyDecWidth

    def _enctab(self, val):
        assert 0 and "not implemented"

    def _dectab(self, seq):
        assert 0 and "not implemented"

    def _encode(self, s):
        key = convertIntToBytes(s[0],self._keyDecWidth)
        payload = bytearray(s[1])
        #print ("{}:{}".format(key,payload))
        enc = []
        for i in range(0,len(key),self._cwDecWidth):
            enc.append( self._enctab( convertBytesToInt(key[i:i+self._cwDecWidth])) )
        for i in range(0,len(payload),self._cwDecWidth):
            enc.append( self._enctab( convertBytesToInt(payload[i:i+self._cwDecWidth])))
        return "".join(enc)

    def _decode(self, s):
        bytes = []
        for i in range(0,len(s),self._cwEncWidth):
            val = self._dectab(s[i:i+self._cwEncWidth])
            val = convertIntToBytes(val,self._cwDecWidth)
            bytes = bytes + val
        key = convertBytesToInt(bytes[0:self._keyDecWidth])
        #bytes = [ chr(b) for b in bytes[self._keyDecWidth:] ]
        return [key,bytes[self._keyDecWidth:]]
        #return [key,bytearray(bytes)]
