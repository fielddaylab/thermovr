using BeauUtil;
using BeauData;
using BeauUtil.Variants;
using System.Collections.Generic;

namespace ThermoVR {
    public class PlayerData : ISerializedObject, ISerializedVersion {

        #region ISerializedObject
        
        // v1: initial
        ushort ISerializedVersion.Version {
            get { return 1; }
        }

        public void Serialize(Serializer serializer) {

        }

        #endregion // ISerializedObject
    }
}